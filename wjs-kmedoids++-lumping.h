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

// Identical hashes for T,U and U,T, but that will never happen since T,U are ordered
struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const pair<T, U> &x) const
  {
    return hash<T>()(x.first) ^ hash<U>()(x.second);
  }
};

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


class MedoidsStateNode{
public:
	MedoidsStateNode();
	MedoidsStateNode(int stateid);
	int stateId;
	bool center = false;
	double minDiv = 1.0;
	int minCenterStateId;
};

MedoidsStateNode::MedoidsStateNode(){
};

MedoidsStateNode::MedoidsStateNode(int stateid){
	stateId = stateid;
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
	vector<int> stateNodeIndices;
	vector<int> stateNodeDanglingIndices;
};

PhysNode::PhysNode(){
};


class StateNetwork{
private:
	void calcEntropyRate();
	double wJSdiv(int stateId1, int stateId2);
	void findCenters(vector<MedoidsStateNode> &medoidsStateNodes);
	void performLumping(vector<MedoidsStateNode> &medoidsStateNodes);

	string inFileName;
	string outFileName;
	std::mt19937 &mtRand;
	int NphysNodes = 0;
	int NstateNodes = 0;
	int Nlinks = 0;
	int Ndanglings = 0;
	int Ncontexts = 0;
	int NphysDanglings = 0;
	int Nclu;
	double totWeight = 0.0;
	double entropyRate = 0.0;
	unordered_map<int,PhysNode> physNodes;
	unordered_map<pair<int,int>,double,pairhash> cachedWJSdiv;
	vector<StateNode> stateNodes;

public:
	StateNetwork(string infilename,string outfilename,int nclu,std::mt19937 &mtrand);
	
	// void lumpDanglings();
	void lumpStateNodes();
	void loadStateNetwork();
	void printStateNetwork();

};

StateNetwork::StateNetwork(string infilename,string outfilename,int nclu,std::mt19937 &mtrand) : mtRand(mtrand){
	Nclu = nclu;
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

	unordered_map<pair<int,int>,double,pairhash>::iterator wJSdiv_it = cachedWJSdiv.find(make_pair(stateIndex1,stateIndex2));
	if(wJSdiv_it != cachedWJSdiv.end())
		return wJSdiv_it->second;

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

		if((links1 != links1end && links2 == links2end) || ((links1 != links1end && links2 != links2end) && (links1->first < links2->first))){
		// If the first state node has a link that the second has not

			double p1 = links1->second/ow1;
			h1 -= p1*log(p1);
			double p12 = pi1*links1->second/ow1;
			h12 -= p12*log(p12);
			links1++;

		}
		else if((links2 != links2end && links1 == links1end) || ((links1 != links1end && links2 != links2end) && (links2->first < links1->first))){
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
			double p2 = links2->second/ow2;
			h2 -= p2*log(p2);
			double p12 = pi1*links1->second/ow1 + pi2*links2->second/ow2;
			h12 -= p12*log(p12);
			links1++;
			links2++;

		}
	}
	double div = (w1+w2)*h12 - w1*h1 - w2*h2;
	cachedWJSdiv[make_pair(stateIndex1,stateIndex2)] = div;
	return div;
}

void StateNetwork::calcEntropyRate(){

	cout << "Calculating entropy rate:" << endl;
	
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

	cout << "-->Entropy rate in bits: " << entropyRate << endl;

}

void StateNetwork::performLumping(vector<MedoidsStateNode> &medoidsStateNodes){

	int NPstateNodes = medoidsStateNodes.size();
	// Update stateNodes to reflect the lumping
	for(int i=Nclu;i<NPstateNodes;i++){

		StateNode &lumpedStateNode = stateNodes[medoidsStateNodes[i].minCenterStateId];
		StateNode &lumpingStateNode = stateNodes[medoidsStateNodes[i].stateId];

		// Add context to lumped state node
		lumpedStateNode.contexts.insert(lumpedStateNode.contexts.begin(),lumpingStateNode.contexts.begin(),lumpingStateNode.contexts.end());
		// Add links to lumped state node
		Nlinks  -= lumpedStateNode.links.size() + lumpingStateNode.links.size(); // To update the global number of links
		for(map<int,double>::iterator link_it = lumpingStateNode.links.begin(); link_it != lumpingStateNode.links.end(); link_it++)
			lumpedStateNode.links[link_it->first] += link_it->second;
		Nlinks += lumpedStateNode.links.size(); // To update the global number of links
		lumpedStateNode.outWeight += lumpingStateNode.outWeight;

		// Update state id of lumping state node to point to lumped state node and make it inactive
		lumpingStateNode.stateId = lumpedStateNode.stateId;
		lumpingStateNode.active = false;

	}

}

void StateNetwork::findCenters(vector<MedoidsStateNode> &medoidsStateNodes){
	// Modifies the order of medoidsStateNodes such thar the fist Nclu will be the centers.
	// Also, all elements will contain the stateId it is closest to.

	int NPstateNodes = medoidsStateNodes.size();
	int Ncenters = 0;
	double sumMinDiv = 1.0*NPstateNodes; // Because minDiv is set to 1.0 for all state nodes

	// Find random state node in physical node as first center
	std::uniform_int_distribution<int> randInt(0,NPstateNodes-1);
	int firstCenterIndex = randInt(mtRand);
	medoidsStateNodes[firstCenterIndex].center = true;
	medoidsStateNodes[firstCenterIndex].minCenterStateId = medoidsStateNodes[firstCenterIndex].stateId;
	sumMinDiv -= medoidsStateNodes[firstCenterIndex].minDiv;
	// Put the center in first non-center position (Ncenters = 0) by swapping elements
	swap(medoidsStateNodes[Ncenters],medoidsStateNodes[firstCenterIndex]);
	Ncenters++;


	// Find Nclu-1 more centers based on the kmedoids++ algorithm
	while(Ncenters < Nclu){
		int lastClusterId = medoidsStateNodes[Ncenters-1].stateId;
		for(int i=Ncenters;i<NPstateNodes;i++){
			double div = wJSdiv(medoidsStateNodes[i].stateId,lastClusterId);
			if(div < medoidsStateNodes[i].minDiv){
				// Found new minimum divergence to center
				sumMinDiv -= medoidsStateNodes[i].minDiv;
				sumMinDiv += div;
				medoidsStateNodes[i].minDiv = div;
				medoidsStateNodes[i].minCenterStateId = lastClusterId;
			}
		}
		// Pick new center proportional to minimum divergence
		std::uniform_int_distribution<double> randDouble(0.0,sumMinDiv);
		double randMinDivSum = randDouble(mtRand);
		double minDivSum = 0.0;
		int newCenterIndex = Ncenters;
		for(int i=Ncenters;i<NPstateNodes;i++){
			minDivSum += medoidsStateNodes[i].minDiv;
			if(minDivSum > randMinDivSum){
				newCenterIndex = i;
				break;
			}
		}
		medoidsStateNodes[newCenterIndex].center = true;
		medoidsStateNodes[newCenterIndex].minCenterStateId = medoidsStateNodes[newCenterIndex].stateId;
		sumMinDiv -= medoidsStateNodes[newCenterIndex].minDiv;
		// Put the center in first non-center position by swapping elements
		swap(medoidsStateNodes[Ncenters],medoidsStateNodes[newCenterIndex]);
		Ncenters++;
	}

	// Check if last center gives minimum divergence for some state nodes
	int lastClusterId = medoidsStateNodes[Ncenters-1].stateId;
		for(int i=Ncenters;i<NPstateNodes;i++){
			double div = wJSdiv(medoidsStateNodes[i].stateId,lastClusterId);
			if(div < medoidsStateNodes[i].minDiv){
				// Found new minimum divergence to center
				sumMinDiv -= medoidsStateNodes[i].minDiv;
				sumMinDiv += div;
				medoidsStateNodes[i].minDiv = div;
				medoidsStateNodes[i].minCenterStateId = lastClusterId;
			}
		}

}

void StateNetwork::lumpStateNodes(){

	cout << "Lumping state nodes in each physical node:" << endl;

	int Nlumpings = 0;
	int Nprocessed = 0;

	for(unordered_map<int,PhysNode>::iterator phys_it = physNodes.begin(); phys_it != physNodes.end(); phys_it++){

		PhysNode &physNode = phys_it->second;
		int NPstateNodes = physNode.stateNodeIndices.size();
		
		if(NPstateNodes > Nclu){

			// Initialize vector with state nodes in physical node with minimum necessary information
			vector<MedoidsStateNode> medoidsStateNodes = vector<MedoidsStateNode>(NPstateNodes);
			for(int i=0;i<NPstateNodes;i++){
				medoidsStateNodes[i].stateId = physNode.stateNodeIndices[i];
			}

			// Initialize vector to store centers and stateId of closest center
			findCenters(medoidsStateNodes);

			// Free cached divergences
			cachedWJSdiv = unordered_map<pair<int,int>,double,pairhash>();

			// Perform the lumping and update stateNodes
			performLumping(medoidsStateNodes);
			Nlumpings += NPstateNodes - Nclu;
			NstateNodes -= NPstateNodes - Nclu;

		}
		Nprocessed++;
		cout << "\r-->Lumped " << Nlumpings << " state nodes in " << Nprocessed << "/" << NphysNodes << " physical nodes.               ";
	}

	cout << endl;

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
				physNodes[physId].stateNodeIndices.push_back(stateId);
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

// void StateNetwork::lumpDanglings(){

// 	unordered_set<int> physDanglings;
// 	int Nlumpings = 0;

// 	cout << "Lumping dangling state nodes:" << endl;

// 	// First loop sets updated stateIds of non-dangling state nodes and state nodes in dangling physical nodes, which are lumped into one state node
// 	int updatedStateId = 0;
// 	for(vector<StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
// 		if(it->outWeight > epsilon){
// 			// Set updated stateIds for non-dangling state nodes
// 			it->stateId = updatedStateId;
// 			updatedStateId++;
// 		}
// 		else{
// 			// Lump all dangling state nodes into one state node in dangling physical nodes, and update the stateIds
// 			int NnonDanglings = physNodes[it->physId].stateNodeIndices.size();
// 			if(NnonDanglings == 0){

// 				// When all state nodes are dangling, lump them to the first dangling state node id of the physical node
// 				physDanglings.insert(it->physId);
// 				// Id of first dangling state node
// 				int lumpedStateIndex = physNodes[it->physId].stateNodeDanglingIndices[0];
// 				if(lumpedStateIndex == it->stateId){
// 					// The first dangling state node in dangling physical node remains
// 					it->stateId = updatedStateId;
// 					updatedStateId++;
// 				}	
// 				else{
// 					// Add context to lumped state node
// 					stateNodes[lumpedStateIndex].contexts.insert(stateNodes[lumpedStateIndex].contexts.begin(),it->contexts.begin(),it->contexts.end());
// 					// Update state id to point to lumped state node with upodated stateId and make it inactive
// 					it->stateId = stateNodes[lumpedStateIndex].stateId;
// 					it->active = false;
// 					// Number of state nodes reduces by 1
// 					NstateNodes--;
// 					Nlumpings++;
// 				}	
// 			}
// 		}
// 	}

// 	// Second loop sets updated stateIds of dangling state nodes in physical nodes with non-dangling state nodes
// 	for(vector<StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){

// 		if(it->outWeight < epsilon){
// 			int NnonDanglings = physNodes[it->physId].stateNodeIndices.size();
// 			if(NnonDanglings > 0){

// 				std::uniform_int_distribution<int> randInt(0,NnonDanglings-1);
// 				// Find random state node
// 				int lumpedStateIndex = physNodes[it->physId].stateNodeIndices[randInt(mtRand)];
// 				// Add context to lumped state node
// 				stateNodes[lumpedStateIndex].contexts.insert(stateNodes[lumpedStateIndex].contexts.begin(),it->contexts.begin(),it->contexts.end());
				
// 				// Update state id to point to lumped state node and make it inactive
// 				it->stateId = stateNodes[lumpedStateIndex].stateId;

// 				it->active = false;
// 				// Number of state nodes reduces by 1
// 				NstateNodes--;
// 				Nlumpings++;

// 			}
// 		}
// 	}

// 	NphysDanglings = physDanglings.size();
// 	cout << "-->Lumped " << Nlumpings << " dangling state nodes." << endl;
// 	cout << "-->Found " << NphysDanglings << " dangling physical nodes. Lumped dangling state nodes into a single dangling state node." << endl;

// 	cout << physNodes[1].stateNodeIndices.size() << endl;
// 	for(int i=0;i<10;i++){
// 		cout << wJSdiv(physNodes[1].stateNodeIndices[0],physNodes[1].stateNodeIndices[i]) << endl;
// 	}
// }