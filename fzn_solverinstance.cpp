/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Guido Tack <guido.tack@monash.edu>
 */

 /* This Source Code Form is subject to the terms of the Mozilla Public
  * License, v. 2.0. If a copy of the MPL was not distributed with this
  * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifdef _WIN32
#define NOMINMAX     // Need this before all (implicit) include's of Windows.h
#endif

#include "minizinc/solvers/fzn_solverinstance.hh"
const auto SolverInstance__ERROR = MiniZinc::SolverInstance::ERROR;  // before windows.h
#include <cstdio>
#include <fstream>



#include <minizinc/timer.hh>
#include <minizinc/prettyprinter.hh>
#include <minizinc/nlprettyprinter.hh>
#include <minizinc/parser.hh>
#include <minizinc/typecheck.hh>
#include <minizinc/builtins.hh>
#include <minizinc/eval_par.hh>

//#include </Users/matthewskelley/libminizinc/lib/nlprettyprinter.cpp>  //Added to use the variableVector


#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#include <tchar.h>
//#include <atlstr.h>
#else
#include <unistd.h>
#include <sys/select.h>
#include <sys/time.h>
#endif
#include <sys/types.h>
#include <signal.h>
#include <thread>
#include <mutex>

namespace MiniZinc {

  class FZN_SolverFactory: public SolverFactory {
    Options _options;
  public:
    SolverInstanceBase* doCreateSI(Env& env) {
      return new FZNSolverInstance(env, _options);
    }
    string getVersion( );
    bool processOption(int& i, int argc, const char** argv);
    void printHelp(std::ostream& os);
    
    
    
  };
// #define __NO_EXPORT_FZN_SOLVERINSTANCE__  // define this to avoid exporting
#ifndef __NO_EXPORT_FZN_SOLVERINSTANCE__
  FZN_SolverFactory fzn_solverfactory;
#endif

  string FZN_SolverFactory::getVersion()
  {
    string v = "NICTA FZN solver plugin, compiled  " __DATE__ "  " __TIME__;
    return v;
  }

  bool FZN_SolverFactory::processOption(int& i, int argc, const char** argv)
  {
    if (string(argv[i])=="--solver") {
      i++;
      if (i==argc) {
        goto error;
      }
      _options.setStringParam(constants().opts.solver.fzn_solver.str(), argv[i]);
    } else if (string(argv[i])=="--use_nl") {
      _options.setBoolParam("use_nl", true);
    }
    return true;
  error:
    return false;
  }
  
  void FZN_SolverFactory::printHelp(ostream& os)
  {
    os
    << "FZN solver plugin options:" << std::endl
    << "--solver         the backend solver"
    << std::endl;
  }

  namespace {
    
//    std::vector<Expression*> genOutputVariables;
//    std::vector<Expression*> genOutputArrays;

#ifdef _WIN32
    mutex mtx;
    void ReadPipePrint(HANDLE g_hCh, ostream& os, ostream* pResult = nullptr) {
      bool done = false;
      while (!done) {
        char buffer[5255];
        DWORD count = 0;
        bool bSuccess = ReadFile(g_hCh, buffer, sizeof(buffer) - 1, &count, NULL);
        if (bSuccess && count > 0) {
          buffer[count] = 0;
          lock_guard<mutex> lck(mtx);
          do {
            if (pResult)
              (*pResult) << buffer;
            os << buffer << flush;
          } while (0);
        }
        else {
          done = true;
        }
      }
    }
#endif

    class FznProcess {
    protected:
      std::string _fzncmd;
      bool _canPipe;
      Model* _flat;
    public:
      FznProcess(const std::string& fzncmd, bool pipe, Model* flat) : _fzncmd(fzncmd), _canPipe(pipe), _flat(flat) {}
      std::string run(bool use_nl) {
#ifdef _WIN32
        std::stringstream result;

        SECURITY_ATTRIBUTES saAttr;
        saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
        saAttr.bInheritHandle = TRUE;
        saAttr.lpSecurityDescriptor = NULL;

        HANDLE g_hChildStd_IN_Rd = NULL;
        HANDLE g_hChildStd_IN_Wr = NULL;
        HANDLE g_hChildStd_OUT_Rd = NULL;
        HANDLE g_hChildStd_OUT_Wr = NULL;
        HANDLE g_hChildStd_ERR_Rd = NULL;
        HANDLE g_hChildStd_ERR_Wr = NULL;

        // Create a pipe for the child process's STDOUT. 
        if (!CreatePipe(&g_hChildStd_OUT_Rd, &g_hChildStd_OUT_Wr, &saAttr, 0))
          std::cerr << "Stdout CreatePipe" << std::endl;
        // Ensure the read handle to the pipe for STDOUT is not inherited.
        if (!SetHandleInformation(g_hChildStd_OUT_Rd, HANDLE_FLAG_INHERIT, 0))
          std::cerr << "Stdout SetHandleInformation" << std::endl;

        // Create a pipe for the child process's STDERR. 
        if (!CreatePipe(&g_hChildStd_ERR_Rd, &g_hChildStd_ERR_Wr, &saAttr, 0))
          std::cerr << "Stderr CreatePipe" << std::endl;
        // Ensure the read handle to the pipe for STDERR is not inherited.
        if (!SetHandleInformation(g_hChildStd_ERR_Rd, HANDLE_FLAG_INHERIT, 0))
          std::cerr << "Stderr SetHandleInformation" << std::endl;

        // Create a pipe for the child process's STDIN
        if (!CreatePipe(&g_hChildStd_IN_Rd, &g_hChildStd_IN_Wr, &saAttr, 0))
          std::cerr << "Stdin CreatePipe" << std::endl;
        // Ensure the write handle to the pipe for STDIN is not inherited. 
        if (!SetHandleInformation(g_hChildStd_IN_Wr, HANDLE_FLAG_INHERIT, 0))
          std::cerr << "Stdin SetHandleInformation" << std::endl;

        std::string fznFile;
        if (!_canPipe) {
          TCHAR szTempFileName[MAX_PATH];
          TCHAR lpTempPathBuffer[MAX_PATH];

          GetTempPath(MAX_PATH, lpTempPathBuffer);
          GetTempFileName(lpTempPathBuffer,
            "tmp_fzn_", 0, szTempFileName);

          fznFile = szTempFileName;
          MoveFile(fznFile.c_str(), (fznFile + ".fzn").c_str());
          fznFile += ".fzn";
          std::ofstream os(fznFile);
          Printer p(os, 0, true);
          p.print(_flat);
        }

        PROCESS_INFORMATION piProcInfo;
        STARTUPINFO siStartInfo;
        BOOL bSuccess = FALSE;

        // Set up members of the PROCESS_INFORMATION structure.
        ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));

        // Set up members of the STARTUPINFO structure. 
        // This structure specifies the STDIN and STDOUT handles for redirection.
        ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));
        siStartInfo.cb = sizeof(STARTUPINFO);
        siStartInfo.hStdError = g_hChildStd_ERR_Wr;
        siStartInfo.hStdOutput = g_hChildStd_OUT_Wr;
        siStartInfo.hStdInput = g_hChildStd_IN_Rd;
        siStartInfo.dwFlags |= STARTF_USESTDHANDLES;

        std::stringstream cmdline;
        cmdline << strdup(_fzncmd.c_str()) << " " << strdup("-v") << " ";
        if (_canPipe) {
          cmdline << strdup("-");
        }
        else {
          cmdline << strdup(fznFile.c_str());
        }

        char* cmdstr = strdup(cmdline.str().c_str());

        bool processStarted = CreateProcess(NULL,
          cmdstr,        // command line 
          NULL,          // process security attributes 
          NULL,          // primary thread security attributes 
          TRUE,          // handles are inherited 
          0,             // creation flags 
          NULL,          // use parent's environment 
          NULL,          // use parent's current directory 
          &siStartInfo,  // STARTUPINFO pointer 
          &piProcInfo);  // receives PROCESS_INFORMATION 

        if (!processStarted) {
//          std::cout<<"!processStarted"<<std::endl;
          std::stringstream ssm;
          ssm << "Error occurred when executing FZN solver with command \"" << cmdstr << "\".";
          throw InternalError(ssm.str());
        }

        CloseHandle(piProcInfo.hProcess);
        CloseHandle(piProcInfo.hThread);
        delete cmdstr;

        if (_canPipe) {
          DWORD dwWritten;
          for (Model::iterator it = _flat->begin(); it != _flat->end(); ++it) {
            std::stringstream ss;
            Item* item = *it;
            ss << *item;
            std::string str = ss.str();
            bSuccess = WriteFile(g_hChildStd_IN_Wr, str.c_str(),
              str.size(), &dwWritten, NULL);
          }
        }

        // Stop ReadFile from blocking
        CloseHandle(g_hChildStd_OUT_Wr);
        CloseHandle(g_hChildStd_ERR_Wr);
        // Just close the child's in pipe here
        CloseHandle(g_hChildStd_IN_Rd);

        // Threaded solution seems simpler than asyncronous pipe reading
        thread thrStdout(ReadPipePrint, g_hChildStd_OUT_Rd, ref(cout), &(result));
        thread thrStderr(ReadPipePrint, g_hChildStd_ERR_Rd, ref(cerr), nullptr);
        thrStdout.join();
        thrStderr.join();

        // Hard timeout: GenerateConsoleCtrlEvent()

        if (!_canPipe) {
          remove(fznFile.c_str());
        }

        return result.str();
      }
      
      
#else
      
      std::vector<Expression*> genOutputVariables;
      std::vector<Expression*> genOutputArrays;
      
        int pipes[3][2];
        pipe(pipes[0]);
        pipe(pipes[1]);
        pipe(pipes[2]);

        std::vector<Expression*> varVect;
        std::string fznFile;
        std::string nlSolFile;
        if (!_canPipe) {
          
          
          
          if(!use_nl) {
            char tmpfile[] = "/tmp/fznfileXXXXXX.fzn";
            mkstemps(tmpfile, 4);
            fznFile = tmpfile;
            std::ofstream os(tmpfile);
            Printer p(os, 0);
            p.print(_flat);
          } else {
            
            /*
             Create the .nl file we will write to
             */
            char tmpfile[] = "/tmp/fznfileXXXXXX.nl";
            mkstemps(tmpfile, 3);
            fznFile = tmpfile;
            std::ofstream os(tmpfile);
            NLPrinter p(os);
            p.print(_flat);
            
            /* 
             Don't need this because the .sol file is created with the same directory as the .nl file
             */
            char nltmpfile[] = "/tmp/nlXXXXXX.sol";
            mkstemps(nltmpfile, 4);
            nlSolFile = nltmpfile;
            
            genOutputVariables = p.outputVariables;
            genOutputArrays = p.outputArrays;
          }
        }

        // Make sure to reap child processes to avoid creating zombies
        signal(SIGCHLD, SIG_IGN);
            
        if (int childPID = fork()) {
          close(pipes[0][0]);
          close(pipes[1][1]);
          close(pipes[2][1]);
          if (_canPipe) {
            for (Model::iterator it = _flat->begin(); it != _flat->end(); ++it) {
              std::stringstream ss;
              Item* item = *it;
              ss << *item;
              std::string str = ss.str();
              write(pipes[0][1], str.c_str(), str.size());
            }
          }
          close(pipes[0][1]);
          std::stringstream result;

          fd_set fdset;
          FD_ZERO(&fdset);

          bool done = false;
          while (!done) {
            FD_SET(pipes[1][0], &fdset);
            FD_SET(pipes[2][0], &fdset);
            if ( 0>=select(FD_SETSIZE, &fdset, NULL, NULL, NULL) )
            {
              kill(childPID, SIGKILL);
              done = true;
            } else {
              for ( int i=1; i<=2; ++i )
                if ( FD_ISSET( pipes[i][0], &fdset ) )
                {
                  char buffer[1000];
                  int count = read(pipes[i][0], buffer, sizeof(buffer) - 1);
                  if (count > 0) {
                    buffer[count] = 0;
                    if ( 1==i ) {
                      result << buffer;
                      cout << buffer << flush;
                    }
                    else
                      cerr << buffer << flush;
                  }
                  else {
                    done = true;
                  }
                }
            }
          }
              
              
          std::string out;
          // if IPOPT : read from .sol file and process it, leaving the result in out

              if(!use_nl) {
                out = result.str();
              } else {
                
                std::ostringstream oss;
                
                string tmp = fznFile.c_str();
                tmp.erase(tmp.end()-2,tmp.end());
                tmp = tmp + "sol";
                ifstream inFile;
                inFile.open(tmp);
                
                string outputCondition;
                string constraints;
                string variables;
                string currentInput;
                vector<string> solVarVector;
                
                /*
                 The following lines read the .sol file and extract the solution
                 */
                
                while (currentInput != "3.11.1:") {
                  inFile>>currentInput;
                }
                
                inFile>>currentInput;
                outputCondition.append(currentInput);
                
                while (currentInput != "Options") {
                  inFile>>currentInput;
                  if (currentInput == "Options") {
                    break;
                  }
                  outputCondition.append(" " + currentInput);
                }
                
                for (int j = 0; j < 6; j++) {
                  inFile>>constraints;
                }
                
                for (int j = 0; j < 2; j++) {
                  inFile>>variables;
                }
                
                for (int j = 0; j <std::atoi(constraints.c_str()); j++) {
                  inFile>>currentInput;
                }
                for (int j = 0; j < std::atoi(variables.c_str()); j++) {
                  inFile>>currentInput;
                  solVarVector.push_back(currentInput);
                }
                /*
                 Prints the output variables
                 */
                for (int i = 0; i < genOutputVariables.size(); i++) {
                  VarDecl* vd = genOutputVariables[i]->cast<VarDecl>();
                  oss<<"\n"<<vd->id()->v()<<" = "<<solVarVector[vd->payload()]<<"; "<<std::endl;
                }
                
                /*
                Prints the output arrays
                */
                for (int i = 0; i < genOutputArrays.size(); i++) {
                  VarDecl* vd = genOutputArrays[i]->cast<VarDecl>();
                  const Expression* e = vd->e();
                  const ArrayLit& al = *e->cast<ArrayLit>();
                  int n = al.dims();
                  oss<<"\n"<<*vd->id()<<" = ";
                  if (n == 1 && al.min(0) == 1) {
                    oss << "[";
                    for (unsigned int i = 0; i < al.v().size(); i++) {
                      oss<<solVarVector[al.v()[i]->cast<Id>()->decl()->payload()];
                      if (i<al.v().size()-1) {
                        oss << ",";
                      }
                    }
                    oss << "]; "<<std::endl;
                  }
                }
                oss << "----------\n" << "==========\n";
                if (outputCondition == "Optimal Solution Found") {
                  oss << "==========\n";
                }
                else {
                  oss << "=====UNSATISFIABLE=====";
                }
                out = oss.str();
              }
          if (!_canPipe) {
            remove(fznFile.c_str());
            if(use_nl)
              remove(nlSolFile.c_str());
          }
          return out;
        }
        else {
          close(STDOUT_FILENO);
          close(STDERR_FILENO);
          close(STDIN_FILENO);
          dup2(pipes[0][0], STDIN_FILENO);
          dup2(pipes[1][1], STDOUT_FILENO);
          dup2(pipes[2][1], STDERR_FILENO);
          close(pipes[0][0]);
          close(pipes[0][1]);
          close(pipes[1][1]);
          close(pipes[1][0]);
          close(pipes[2][1]);
          close(pipes[2][0]);

          std::vector<char*> cmd_line;
          
          if(!use_nl) {
            cmd_line.push_back(strdup(_fzncmd.c_str()));
            cmd_line.push_back(strdup(_canPipe ? "-" : fznFile.c_str()));
            cmd_line.push_back(strdup("-v"));
          } else {
            cmd_line.push_back(strdup("ipopt"));
            cmd_line.push_back(strdup("-s"));
            cmd_line.push_back(strdup(fznFile.c_str()));

          }

          char** argv = new char*[cmd_line.size() + 1];
          for (unsigned int i = 0; i < cmd_line.size(); i++)
            argv[i] = cmd_line[i];
          argv[cmd_line.size()] = 0;
//          std::cerr << "CMD: \"" << argv[0] << " " << argv[1] << " " << argv[2] << "\".";

          int status = execvp(argv[0], argv);
          if (status == -1) {
            std::stringstream ssm;
            ssm << "Error occurred when executing FZN solver with command \"" << argv[0] << " " << argv[1] << " " << argv[2] << "\".";
            throw InternalError(ssm.str());
          }
        }
        assert(false);
    }
#endif
    };
  }

  FZNSolverInstance::FZNSolverInstance(Env& env, const Options& options)
    : SolverInstanceImpl<FZNSolver>(env, options), _fzn(env.flat()), _ozn(env.output()), hadSolution(false) {}

  FZNSolverInstance::~FZNSolverInstance(void) {}

  namespace {
    ArrayLit* b_arrayXd(Env& env, ASTExprVec<Expression> args, int d) {
      GCLock lock;
      ArrayLit* al = eval_array_lit(env.envi(), args[d]);
      std::vector<std::pair<int, int> > dims(d);
      unsigned int dim1d = 1;
      for (int i = 0; i < d; i++) {
        IntSetVal* di = eval_intset(env.envi(), args[i]);
        if (di->size() == 0) {
          dims[i] = std::pair<int, int>(1, 0);
          dim1d = 0;
        }
        else if (di->size() != 1) {
          throw EvalError(env.envi(), args[i]->loc(), "arrayXd only defined for ranges");
        }
        else {
          dims[i] = std::pair<int, int>(static_cast<int>(di->min(0).toInt()),
            static_cast<int>(di->max(0).toInt()));
          dim1d *= dims[i].second - dims[i].first + 1;
        }
      }
      if (dim1d != al->v().size())
        throw EvalError(env.envi(), al->loc(), "mismatch in array dimensions");
      ArrayLit* ret = new ArrayLit(al->loc(), al->v(), dims);
      Type t = al->type();
      t.dim(d);
      ret->type(t);
      ret->flat(al->flat());
      return ret;
    }
  }

  SolverInstance::Status
    FZNSolverInstance::solve(void) {
    std::vector<std::string> includePaths;
    std::string fzn_solver = _options.getStringParam(constants().opts.solver.fzn_solver.str(), "flatzinc");
    bool use_nl = _options.getBoolParam("use_nl", false);
    if (_options.getBoolParam(constants().opts.verbose.str(), false)) {
      std::cerr << "Using FZN solver " << fzn_solver << " for solving." << std::endl;
    }
    FznProcess proc(fzn_solver, false, _fzn);
    std::string r = proc.run(use_nl);
    std::stringstream result;
    result << r;
    std::string solution;

    typedef std::pair<VarDecl*, Expression*> DE;
    ASTStringMap<DE>::t declmap;
    for (unsigned int i = 0; i < _ozn->size(); i++) {
      if (VarDeclI* vdi = (*_ozn)[i]->dyn_cast<VarDeclI>()) {
        declmap.insert(std::make_pair(vdi->e()->id()->str(), DE(vdi->e(), vdi->e()->e())));
      }
    }

    hadSolution = false;
    while (result.good()) {
      std::string line;
      getline(result, line);

      if (beginswith(line, "----------")) {
        if (hadSolution) {
          for (ASTStringMap<DE>::t::iterator it = declmap.begin(); it != declmap.end(); ++it) {
            it->second.first->e(it->second.second);
          }
        }
        Model* sm = parseFromString(solution, "solution.szn", includePaths, true, false, false, std::cerr);
        if (sm) {
          for (Model::iterator it = sm->begin(); it != sm->end(); ++it) {
            if (AssignI* ai = (*it)->dyn_cast<AssignI>()) {
              ASTStringMap<DE>::t::iterator it = declmap.find(ai->id());
              if (it == declmap.end()) {
                std::cerr << "Error: unexpected identifier " << ai->id() << " in output\n";
                exit(EXIT_FAILURE);
              }
              if (Call* c = ai->e()->dyn_cast<Call>()) {
                // This is an arrayXd call, make sure we get the right builtin
                assert(c->args()[c->args().size() - 1]->isa<ArrayLit>());
                for (unsigned int i = 0; i < c->args().size(); i++)
                  c->args()[i]->type(Type::parsetint());
                c->args()[c->args().size() - 1]->type(it->second.first->type());
                ArrayLit* al = b_arrayXd(_env, c->args(), c->args().size() - 1);
                it->second.first->e(al);
              }
              else {
                it->second.first->e(ai->e());
              }
            }
          }
          delete sm;
          hadSolution = true;
        } else {
          std::cerr << "\n\n\nError: solver output malformed; DUMPING: ---------------------\n\n\n";
          cerr << result.str() << endl;
          std::cerr << "\n\nError: solver output malformed (END DUMPING) ---------------------" << endl;
          exit(EXIT_FAILURE);
        }
      }
      else if (beginswith(line, "==========")) {
        return hadSolution ? SolverInstance::OPT : SolverInstance::UNSAT;
      }
      else if (beginswith(line, "=====UNSATISFIABLE=====")) {
        return SolverInstance::UNSAT;
      }
      else if (beginswith(line, "=====UNBOUNDED=====")) {
        return SolverInstance::UNBND;
      }
      else if (beginswith(line, "=====UNSATorUNBOUNDED=====")) {
        return SolverInstance::UNSATorUNBND;
      }
      else if (beginswith(line, "=====ERROR=====")) {
        return SolverInstance__ERROR;
      }
      else if (beginswith(line, "=====UNKNOWN=====")) {
        return SolverInstance::UNKNOWN;
      }
      else {
        solution += line;
      }
    }

    return hadSolution ? SolverInstance::SAT : SolverInstance::UNSAT;
  }

  void FZNSolverInstance::printSolution(ostream& os) {
    assert(hadSolution);
    _env.evalOutput(os);
  }

  void
    FZNSolverInstance::processFlatZinc(void) {}

  void
    FZNSolverInstance::resetSolver(void) {}

  Expression*
    FZNSolverInstance::getSolutionValue(Id* id) {
    assert(false);
    return NULL;
  }
}
