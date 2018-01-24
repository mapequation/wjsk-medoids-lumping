#include "wjs-kmedoids++-lumping.h"

using namespace std;
using std::cout;
using std::cin;
using std::endl;

unsigned stou(char *s){
  return strtoul(s,(char **)NULL,10);
}

  // Call: trade <seed> <Ntries>
int main(int argc,char *argv[]){

  cout << "Version: July 17, 2016." << endl;
  cout << "Command: ";
  cout << argv[0];
  for(int i=1;i<argc; i++)
    cout << " " << argv[i];
  cout << endl;

  // Parse command input
  const string CALL_SYNTAX = "Call: ./dangling-lumping [-s <seed>] [-N <number of attempts>] [-k <max number of clusters>] [-e <max entropy in lumped state node>] [--max-order <max Markov order>] [-d <number of clusters in each division (>= 2)>] [--fast] [--batchoutput] [--state-containers state_assignments.txt] [--context-containers lumped_state_network.net] input_state_network.net output_state_network.net [output_state_container.txt]\n";
  if( argc == 1 ){
    cout << CALL_SYNTAX;
    exit(-1);
  }
  unsigned int seed = 1234;

  string inFileName;
  string outFileName;
  string containerOutFileName = "";
  string containerInFileName = "";
  string contextInFileName = "";

  int argNr = 1;
  unsigned int NfinalClu = 100;
  unsigned int NsplitClu = 2;
  double maxEntropy = -1.0;
  int Nattempts = 1;
  int order = -1;
  bool batchOutput = false;
  bool fast = false;
  while(argNr < argc){
    if(to_string(argv[argNr]) == "-h"){
      cout << CALL_SYNTAX;
      exit(-1);
    }
    else if(to_string(argv[argNr]) == "-s"){
      argNr++;
      seed = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "--batchoutput"){
      batchOutput = true;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "--fast"){
      fast = true;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-N"){
      argNr++;
      Nattempts = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "--max-order"){
      argNr++;
      order = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "--context-containers"){
      argNr++;
      contextInFileName = string(argv[argNr]);
      argNr++;
      containerInFileName = "";
    }   
    else if(to_string(argv[argNr]) == "--state-containers"){
      argNr++;
      containerInFileName = string(argv[argNr]);
      argNr++;
      contextInFileName = "";
    }  
    else if(to_string(argv[argNr]) == "-k"){
      argNr++;
      NfinalClu = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-e"){
      argNr++;
      maxEntropy = atof(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-d"){
      argNr++;
      NsplitClu = atoi(argv[argNr]);
      if(NsplitClu < 2){
        cout << "Command error: -d must be integer larger or equal to 2." << endl;
        cout << CALL_SYNTAX;
        exit(-1);
      }
      argNr++;
    }
    else{

      if(argv[argNr][0] == '-'){
        cout << "Unknown command: " << to_string(argv[argNr]) << endl;
        cout << CALL_SYNTAX;
        exit(-1);
      }

      inFileName = string(argv[argNr]);
      argNr++;
      outFileName = string(argv[argNr]);
      argNr++;
      if(argNr < argc){
        containerOutFileName = string(argv[argNr]);
        argNr++;
      }
    }

  }
  
  cout << "Setup:" << endl;
  cout << "-->Using seed: " << seed << endl;
  cout << "-->Will lump state nodes into max number of clusters per container node: " << NfinalClu << endl;
  if(maxEntropy > 0.0)
    cout << "-->Will lump state nodes into max entropy per state node: " << maxEntropy << endl;
  cout << "-->Will iteratively divide worst cluster into number of clusters: " << NsplitClu << endl;
  cout << "-->Will make number of attempts: " << Nattempts << endl;
  if(order > 1)
    cout << "-->Will first perform context lumping to Markov order: " << order << endl;
  else
    cout << "-->Will not perform any context lumping." << endl;
  if(fast)
    cout << "-->Will use medoid center to approximate cluster entropy rate during assignment." << endl;
  else
    cout << "-->Will use aggregate cluster to calculate entropy rate during assignment." << endl;
  if(contextInFileName != "")
    cout << "-->Will lump higher-order states in context containers of lumped state network file: " << contextInFileName << endl;
  if(containerInFileName != "")
    cout << "-->Will lump states in containers of state-to-container assignment file: " << containerInFileName << endl;
  // if(tune)
  //   cout << "-->Will tune medoids for bestter accuracy." << endl;
  // else
  //   cout << "-->Will not tune medoids for better accuracy." << endl;
  cout << "-->Will read state network from file: " << inFileName << endl;
  cout << "-->Will write processed state network to file: " << outFileName << endl;
  if(containerOutFileName != "")
    cout << "-->Will write state node container assignments to: " << containerOutFileName << endl;

  StateNetwork statenetwork(inFileName,outFileName,contextInFileName,containerInFileName,containerOutFileName,NfinalClu,maxEntropy,NsplitClu,Nattempts,order,fast,batchOutput,seed);

  int NprocessedBatches = 0;
  while(statenetwork.loadStateNetworkBatch()){ // NprocessedBatches < 5 && 
    statenetwork.lumpStateNodes();
    NprocessedBatches++;
    if(statenetwork.keepReading || statenetwork.Nbatches > 1){
      statenetwork.printStateNetworkBatch();
      statenetwork.concludeBatch();
    }
    else{
      statenetwork.printStateNetwork();
      statenetwork.printStateNodeContainer();
      break;
    }
  }

  if(statenetwork.Nbatches > 1){
    statenetwork.compileBatches();
    statenetwork.printStateNodeContainer();
  }

}