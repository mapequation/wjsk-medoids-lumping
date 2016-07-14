#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>
#include <functional>
#include <map>
#include <unordered_map>
#include <unordered_set>
using namespace std;
const double epsilon = 1e-15;

unsigned stou(char *s);

template <class T>
inline std::string to_string (const T& t){
	std::stringstream ss;
	ss << t;
	return ss.str();
}

vector<string> tokenize(const string& str,string& delimiters){

	vector<string> tokens;

  // skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while(string::npos != pos || string::npos != lastPos){

    // found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));

    // skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);

    // find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}

	return tokens;

}


class StateNode{
public:
	StateNode();
	StateNode(int stateid, int physid, double outweight);
	int stateId;
	int physId;
	double outWeight;
	bool active = true;
	map<int,double> links;
	vector<string> contexts;
};

StateNode::StateNode(){
};

StateNode::StateNode(int stateid, int physid, double outweight){
	stateId = stateid;
	physId = physid;
	outWeight = outweight;
}

class PhysNode{
public:
	PhysNode();
	vector<int> stateNodeNonDanglingIndices;
	vector<int> stateNodeDanglingIndices;
};

PhysNode::PhysNode(){
};


class StateNetwork{
private:
	void calcEntropyRate();
	double wJSdiv(int stateId1, int stateId2);
	string inFileName;
	string outFileName;
	std::mt19937 &mtRand;
	int NphysNodes = 0;
	int NstateNodes = 0;
	int Nlinks = 0;
	int Ndanglings = 0;
	int Ncontexts = 0;
	int NphysDanglings = 0;
	double totWeight = 0.0;
	double entropyRate = 0.0;
	unordered_map<int,PhysNode> physNodes;
	vector<StateNode> stateNodes;

public:
	StateNetwork(string infilename,string outfilename,std::mt19937 &mtrand);
	
	void lumpDanglings();
	void loadStateNetwork();
	void printStateNetwork();

};

StateNetwork::StateNetwork(string infilename,string outfilename,std::mt19937 &mtrand) : mtRand(mtrand){
	inFileName = infilename;
	outFileName = outfilename;
	mtRand = mtrand;
}

double StateNetwork::wJSdiv(int stateIndex1, int stateIndex2){

	double h1 = 0.0; // The entropy rate of the first state node
	double h2 = 0.0; // The entropy rate of the second state node
	double h12 = 0.0; // The entropy rate of the lumped state node

	if(stateIndex1 > stateIndex2) // Swap to make stateIndex1 lowest 
		swap(stateIndex1,stateIndex2);

	// The out-link weights of the state nodes
	double ow1 = stateNodes[stateIndex1].outWeight;
	double ow2 = stateNodes[stateIndex2].outWeight;
	// Normalized weights over entire network
	double w1 = ow1/totWeight;
	double w2 = ow2/totWeight;
	// Normalized weights over state nodes 1 and 2
	double pi1 = w1 / (w1 + w2);
	double pi2 = w2 / (w1 + w2);

	if(ow1 < epsilon || ow2 < epsilon){ // If one or both state nodes are dangling
		return 0.0;
	}

	map<int,double>::iterator links1 = stateNodes[stateIndex1].links.begin();
	map<int,double>::iterator links2 = stateNodes[stateIndex2].links.begin();
	map<int,double>::iterator links1end = stateNodes[stateIndex1].links.end();
	map<int,double>::iterator links2end = stateNodes[stateIndex2].links.end();
	
	while(links1 != links1end || links2 != links2end){

		if((links1 != links1end && links2 == links2end) || ((links2 != links2end) && (links1->first < links2->first))){
		// If the first state node has a link that the second has not

			double p1 = links1->second/ow1;
			h1 -= p1*log(p1);
			double p12 = pi1*links1->second/ow1;
			h12 -= p12*log(p12);
			links1++;

		}
		else if((links2 != links2end && links1 == links1end) || ((links1 != links1end) && (links2->first < links1->first))){
		// If the second state node has a link that the second has not

			double p2 = links2->second/ow2;
			h2 -= p2*log(p2);
			double p12 = pi2*links2->second/ow2;
			h12 -= p12*log(p12);
			links2++;

		}
		else{ // If both state nodes have the link

			double p1 = links1->second/ow1;
			h1 -= p1*log(p1);
			links1++;
			double p2 = links2->second/ow2;
			h2 -= p2*log(p2);
			links2++;
			double p12 = pi1*links1->second/ow1 + pi2*links2->second/ow2;
			h12 -= p12*log(p12);

		}
	}

	return (w1+w2)*h12 - w1*h1 - w2*h2;
}

void StateNetwork::calcEntropyRate(){
	
	entropyRate = 0.0;

	for(vector<StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		if(it->active){
			double h = 0.0;
			for(map<int,double>::iterator it_link = it->links.begin(); it_link != it->links.end(); it_link++){
				double p = it_link->second/it->outWeight;
				h -= p*log(p);
			}
			entropyRate += it->outWeight/totWeight*h/log(2.0);
		}
	}

}

void StateNetwork::lumpDanglings(){

	unordered_set<int> physDanglings;
	int Nlumpings = 0;

	cout << "Lumping dangling state nodes:" << endl;

	// First loop sets updated stateIds of non-dangling state nodes and state nodes in dangling physical nodes, which are lumped into one state node
	int updatedStateId = 0;
	for(vector<StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		if(it->outWeight > epsilon){
			// Set updated stateIds for non-dangling state nodes
			it->stateId = updatedStateId;
			updatedStateId++;
		}
		else{
			// Lump all dangling state nodes into one state node in dangling physical nodes, and update the stateIds
			int NnonDanglings = physNodes[it->physId].stateNodeNonDanglingIndices.size();
			if(NnonDanglings == 0){

				// When all state nodes are dangling, lump them to the first dangling state node id of the physical node
				physDanglings.insert(it->physId);
				// Id of first dangling state node
				int lumpedStateIndex = physNodes[it->physId].stateNodeDanglingIndices[0];
				if(lumpedStateIndex == it->stateId){
					// The first dangling state node in dangling physical node remains
					it->stateId = updatedStateId;
					updatedStateId++;
				}	
				else{
					// Add context to lumped state node
					stateNodes[lumpedStateIndex].contexts.insert(stateNodes[lumpedStateIndex].contexts.begin(),it->contexts.begin(),it->contexts.end());
					// Update state id to point to lumped state node with upodated stateId and make it inactive
					it->stateId = stateNodes[lumpedStateIndex].stateId;
					it->active = false;
					// Number of state nodes reduces by 1
					NstateNodes--;
					Nlumpings++;
				}	
			}
		}
	}

	// Second loop sets updated stateIds of dangling state nodes in physical nodes with non-dangling state nodes
	for(vector<StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){

		if(it->outWeight < epsilon){
			int NnonDanglings = physNodes[it->physId].stateNodeNonDanglingIndices.size();
			if(NnonDanglings > 0){

				std::uniform_int_distribution<int> randInt(0,NnonDanglings-1);
				// Find random state node
				int lumpedStateIndex = physNodes[it->physId].stateNodeNonDanglingIndices[randInt(mtRand)];
				// Add context to lumped state node
				stateNodes[lumpedStateIndex].contexts.insert(stateNodes[lumpedStateIndex].contexts.begin(),it->contexts.begin(),it->contexts.end());
				
				// Update state id to point to lumped state node and make it inactive
				it->stateId = stateNodes[lumpedStateIndex].stateId;

				it->active = false;
				// Number of state nodes reduces by 1
				NstateNodes--;
				Nlumpings++;

			}
		}
	}

	NphysDanglings = physDanglings.size();
	cout << "-->Lumped " << Nlumpings << " dangling state nodes." << endl;
	cout << "-->Found " << NphysDanglings << " dangling physical nodes. Lumped dangling state nodes into a single dangling state node." << endl;

	cout << physNodes[1].stateNodeNonDanglingIndices.size() << endl;
	for(int i=0;i<10;i++){
		cout << wJSdiv(physNodes[1].stateNodeNonDanglingIndices[0],physNodes[1].stateNodeNonDanglingIndices[i]) << endl;
	}
}

void StateNetwork::loadStateNetwork(){

	string line;
	string buf;
	istringstream ss;

  // ************************* Read state network ************************* //
	ifstream ifs(inFileName.c_str());
	if(!ifs){
		cout << "failed to open \"" << inFileName << "\" exiting..." << endl;
		exit(-1);
	}
	else{
		cout << "Reading " << inFileName << ":" << endl;
	}

	// Skip header lines starting with #
	while(getline(ifs,line)){
		if(line[0] != '#')
			break;
	}

	// ************************* Read states ************************* //
	ss.clear();
	ss.str(line);
	ss >> buf;
	if(buf != "*States"){
		cout << "Expected *States but read " << buf << ", exiting..." << endl;
		exit(-1);
	}
	ss >> buf;
	NstateNodes = atoi(buf.c_str());
	cout << "-->Reading " << NstateNodes  << " state nodes..." << flush;
	stateNodes = vector<StateNode>(NstateNodes);
	for(int i=0;i<NstateNodes;i++){
		getline(ifs,line);
		if(line[0] != '#'){
			ss.clear();
			ss.str(line);
			ss >> buf;
			int stateId = atoi(buf.c_str());
			ss >> buf;
			int physId = atoi(buf.c_str());
	  	ss >> buf;
	  	double outWeight = atof(buf.c_str());
	  	totWeight += outWeight;
			if(outWeight > epsilon)
				physNodes[physId].stateNodeNonDanglingIndices.push_back(stateId);
			else{
				physNodes[physId].stateNodeDanglingIndices.push_back(stateId);
				Ndanglings++;
			}
			stateNodes[stateId] = StateNode(stateId,physId,outWeight);
		}
		else{
			// One extra step for each # comment.
			i--;
		}
	}
	NphysNodes = physNodes.size();
	cout << "found " << Ndanglings << " dangling state nodes in " << NphysNodes << " physical nodes, done!" << endl;

	// ************************* Read arcs ************************* //
	getline(ifs,line);
	ss.clear();
	ss.str(line);
	ss >> buf;
	if(buf != "*Arcs" && buf != "*Links"){
		cout << "Expected *Arcs or *Links but read " << buf << ", exiting..." << endl;
		exit(-1);
	}
	ss >> buf;
	Nlinks = atoi(buf.c_str());
	cout << "-->Reading " << Nlinks  << " state links..." << flush;


	for(int i=0;i<Nlinks;i++){
		getline(ifs,line);
		if(line[0] != '#'){
			ss.clear();
			ss.str(line);
			ss >> buf;
			int source = atoi(buf.c_str());
			ss >> buf;
			int target = atoi(buf.c_str());
			ss >> buf;
			double linkWeight = atof(buf.c_str());
			stateNodes[source].links[target] += linkWeight;
 		}
 		else{
 			// One extra step for each # comment.
 			i--;
 		}
	}
 	cout << "done!" << endl;

	// ************************* Read contexts ************************* //
	getline(ifs,line);
	ss.clear();
	ss.str(line);
	ss >> buf;
	if(buf != "*Contexts" && buf != "*MemoryNodes"){
		cout << "Expected *Contexts or *MemoryNodes but read " << buf << ", exiting..." << endl;
		exit(-1);
	}
	cout << "-->Reading state node contexts..." << flush;
	while(getline(ifs,line)){
		if(line[0] != '#'){
			ss.clear();
			ss.str(line);
			ss >> buf;
			int stateNodeId = atoi(buf.c_str());
			string context = line.substr(buf.length()+1);
			stateNodes[stateNodeId].contexts.push_back(context);
			Ncontexts++;
 		}
	}
 	cout << "found " << Ncontexts << ", done!" << endl;

}

void StateNetwork::printStateNetwork(){

  calcEntropyRate();

	cout << "Writing to " << outFileName << ":" << endl;
  ofstream ofs(outFileName);
  cout << "-->Writing header comments..." << flush;
  ofs << "# Number of physical nodes: " << NphysNodes << "\n";
  ofs << "# Number of state nodes: " << NstateNodes << "\n";
  ofs << "# Number of dangling physical (and state) nodes: " << NphysDanglings << "\n";
  ofs << "# Number of links: " << Nlinks << "\n";
  ofs << "# Number of contexts: " << Ncontexts << "\n";
  ofs << "# Total weight: " << totWeight << "\n";
  ofs << "# Entropy rate: " << entropyRate << "\n";
	cout << "done!" << endl;

	cout << "-->Writing " << NstateNodes << " state nodes..." << flush;
	ofs << "*States " << NstateNodes << "\n";
	ofs << "#stateId ==> (physicalId, outWeight)\n";
	int index = 0;
	for(vector<StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		if(it->active){
			// The state node has not been lumped to another node (but other nodes may have been lumped to it)
			ofs << it->stateId << " " << it->physId << " " << it->outWeight << "\n";
		}
		index++;
	}
	cout << "done!" << endl;

	cout << "-->Writing " << Nlinks << " links..." << flush;
	ofs << "*Links " << Nlinks << "\n";
	ofs << "#(source target) ==> weight\n";
	for(vector<StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		if(it->active){
			// The state node has not been lumped to another node (but other nodes may have been lumped to it)
			for(map<int,double>::iterator it_link = it->links.begin(); it_link != it->links.end(); it_link++){
				ofs << it->stateId << " " << stateNodes[it_link->first].stateId << " " << it_link->second << "\n";
			}
		}
	}
	cout << "done!" << endl;

	cout << "-->Writing " << Ncontexts << " contexts..." << flush;
	ofs << "*Contexts \n";
	ofs << "#stateId <== (physicalId priorId [history...])\n";
	index = 0;
	for(vector<StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		if(it->active){
		// The state node has not been lumped to another node (but other nodes may have been lumped to it)
			for(vector<string>::iterator it_context = it->contexts.begin(); it_context != it->contexts.end(); it_context++){
				ofs << it->stateId << " " << (*it_context) << "\n";
			}
		}
		index++;
	}
	cout << "done!" << endl;

}