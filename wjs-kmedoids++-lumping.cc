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

  cout << "Version: July 10, 2016." << endl;
  cout << "Command: ";
  cout << argv[0];
  for(int i=1;i<argc; i++)
    cout << " " << argv[i];
  cout << endl;

  // Parse command input
  const string CALL_SYNTAX = "Call: ./dangling-lumping [-s <seed>] -k <number of clusters> -batchoutput input_state_network.net output_state_network.net\n";
  if( argc == 1 ){
    cout << CALL_SYNTAX;
    exit(-1);
  }

  unsigned int seed = 1234;

  string inFileName;
  string outFileName;

  int argNr = 1;
  int Nclu = 100;
  double batchOutput = false;
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
    else if(to_string(argv[argNr]) == "-batchoutput"){
      batchOutput = true;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-k"){
      argNr++;
      Nclu = atoi(argv[argNr]);
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
    }

  }

  cout << "Setup:" << endl;
  cout << "-->Using seed: " << seed << endl;
  cout << "-->Will lump state nodes into (at most) number of clusters per physical node: " << Nclu << endl;
  cout << "-->Will read state network from file: " << inFileName << endl;
  cout << "-->Will write processed state network to file: " << outFileName << endl;

  mt19937 mtRand(seed);

  StateNetwork statenetwork(inFileName,outFileName,Nclu,batchOutput,mtRand);

  while(statenetwork.loadStateNetworkBatch()){
    statenetwork.lumpStateNodes();
    if(statenetwork.keepReading || statenetwork.Nbatches > 1){
      statenetwork.printStateNetworkBatch();
      statenetwork.concludeBatch();
    }
    else{
      statenetwork.printStateNetwork();
      break;
    }
  }

  if(statenetwork.Nbatches > 1)
    statenetwork.compileBatches();

}