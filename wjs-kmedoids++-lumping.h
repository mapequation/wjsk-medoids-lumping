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
#include <limits>
const double epsilon = 1.0e-15;
const double bignum = 1.0;
const double threshold = 1.0e-10;

// ofstream with higher precision to avoid truncation errors
struct my_ofstream : std::ofstream {
  explicit my_ofstream(std::streamsize prec = 15)
  {
    this->precision(prec);
  }
};

enum WriteMode { STATENODES, LINKS, CONTEXTS };

template <class T>
inline std::string to_string (const T& t){
	std::stringstream ss;
	ss << t;
	return ss.str();
}

struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const pair<T, U> &x) const
  {
    return x.first*31 + x.second;
  }
};

// // Identical hashes for T,U and U,T, but that will never happen since T,U are ordered
// struct pairhash {
// public:
//   template <typename T, typename U>
//   std::size_t operator()(const pair<T, U> &x) const
//   {
//     return hash<T>()(x.first) ^ hash<U>()(x.second);
//   }
// };

class LocalStateNode{
public:
	LocalStateNode();
	LocalStateNode(int stateid);
	int stateId;
	double minDiv = bignum;
	int minCenterStateId;
};

LocalStateNode::LocalStateNode(){
};

LocalStateNode::LocalStateNode(int stateid){
	stateId = stateid;
}

class StateNode{
public:
	StateNode();
	StateNode(int stateid, int physid, double outweight);
	int stateId;
	int updatedStateId;
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
	double calcEntropyRate();
	double wJSdiv(int stateId1, int stateId2);
	double findCenters(vector<LocalStateNode> &localStateNodes);
	double updateCenters(vector<LocalStateNode> &localStateNodes);
	void performLumping(vector<LocalStateNode> &localStateNodes);
	bool readLines(string &line,vector<string> &lines);
	void writeLines(ifstream &ifs_tmp, ofstream &ofs, WriteMode &writeMode, string &line,int &batchNr);
	void writeLines(ifstream &ifs_tmp, ofstream &ofs, WriteMode &writeMode, string &line);

	// For all batches
	string inFileName;
	string outFileName;
	string tmpOutFileName;
	string tmpOutFileNameStates;
	string tmpOutFileNameLinks;
	string tmpOutFileNameContexts;
	bool batchOutput = false;
	mt19937 &mtRand;
	ifstream ifs;
  string line = "First line";
  double totWeight = 0.0;
  double accumWeight = 0.0;
  int updatedStateId = 0;
	double entropyRate = 0.0;
	unordered_map<int,int> completeStateNodeIdMapping;
  int totNphysNodes = 0;
	int totNstateNodes = 0;
	int totNlinks = 0;
	int totNdanglings = 0;
	int totNcontexts = 0;
	int totNphysDanglings = 0;

	// For each batch
	double weight = 0.0;
	int NphysNodes = 0;
	int NstateNodes = 0;
	int Nlinks = 0;
	int Ndanglings = 0;
	int Ncontexts = 0;
	int NphysDanglings = 0;
	int Nclu;
	unordered_map<pair<int,int>,double,pairhash> cachedWJSdiv;
	unordered_map<int,int> stateNodeIdMapping;
	unordered_map<int,PhysNode> physNodes;
	unordered_map<int,StateNode> stateNodes;

public:
	StateNetwork(string infilename,string outfilename,int nclu,bool batchoutput,std::mt19937 &mtrand);
	
	void lumpStateNodes();
	bool loadStateNetworkBatch();
	void printStateNetworkBatch();
	void printStateNetwork();
	void concludeBatch();
	void compileBatches();

	bool keepReading = true;
  int Nbatches = 0;

};

StateNetwork::StateNetwork(string infilename,string outfilename,int nclu,bool batchoutput,std::mt19937 &mtrand) : mtRand(mtrand){
	Nclu = nclu;
	inFileName = infilename;
	outFileName = outfilename;
	tmpOutFileName = string(outFileName).append("_tmp");
	tmpOutFileNameStates = string(outFileName).append("_tmpstates");
	tmpOutFileNameLinks = string(outFileName).append("_tmplinks");
	tmpOutFileNameContexts = string(outFileName).append("_tmpcontexts");
	mtRand = mtrand;
	batchOutput = batchoutput;

	// Open state network
	ifs.open(inFileName.c_str());
	if(!ifs){
		cout << "failed to open \"" << inFileName << "\" exiting..." << endl;
		exit(-1);
	}
}

double StateNetwork::wJSdiv(int stateIndex1, int stateIndex2){

	double h1 = 0.0; // The entropy rate of the first state node
	double h2 = 0.0; // The entropy rate of the second state node
	double h12 = 0.0; // The entropy rate of the lumped state node

	if(stateIndex1 == stateIndex2)
		return 0.0;

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

	if(div < epsilon)
		div = epsilon;

	cachedWJSdiv[make_pair(stateIndex1,stateIndex2)] = div;


	return div;
}

double StateNetwork::calcEntropyRate(){

	double h = 0.0;

	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(stateNode.active){
			double H = 0.0;
			for(map<int,double>::iterator it_link = stateNode.links.begin(); it_link != stateNode.links.end(); it_link++){
				double p = it_link->second/stateNode.outWeight;
				H -= p*log(p);
			}
			h += stateNode.outWeight*H/totWeight/log(2.0);
		}
	}

	return h;

}



double StateNetwork::updateCenters(vector<LocalStateNode> &localStateNodes){

	int NPstateNodes = localStateNodes.size();
	// To treat clusters one at a time
	unordered_map<int,vector<int> > medoids;
	medoids.reserve(Nclu);
	for(int i=0;i<NPstateNodes;i++)
		medoids[localStateNodes[i].minCenterStateId].push_back(i);
	
	// Loop over all medoids to find new centers
	for(unordered_map<int,vector<int> >::iterator it = medoids.begin(); it != medoids.end(); it++){
		vector<int> &medoid = it->second;
		int NstatesInMedoid = medoid.size();
		double minDivSumInMedoid = bignum*NstatesInMedoid;
		int minDivSumInMedoidIndex = medoid[0];
		// Find total divergence for each state node in medoid
		for(int j=0;j<NstatesInMedoid;j++){
			double divSumInMedoid = 0.0;
			for(int k=0;k<NstatesInMedoid;k++){
				if(k != j){
					divSumInMedoid += wJSdiv(localStateNodes[medoid[j]].stateId,localStateNodes[medoid[k]].stateId);
				}
			}
			if(divSumInMedoid < minDivSumInMedoid){
				minDivSumInMedoid = divSumInMedoid;
				minDivSumInMedoidIndex = j;
			}
		}
		// Update localStateNodes to have all minCenterStateId point to the global Id/Index and to keep centers within the first Nclu elements.
		for(int j=0;j<NstatesInMedoid;j++){
			localStateNodes[medoid[j]].minCenterStateId = localStateNodes[medoid[minDivSumInMedoidIndex]].stateId;
		}
		swap(localStateNodes[medoid[0]],localStateNodes[medoid[minDivSumInMedoidIndex]]); // Swap such that the first Nclu state nodes in localStateNodes are centers
		swap(medoid[0],medoid[minDivSumInMedoidIndex]); // Swap such that the first index in medoid is center
	}

	double sumMinDiv = bignum*(NPstateNodes-Nclu);
	// Find closest center for all state nodes in physical node
	for(unordered_map<int,vector<int> >::iterator it = medoids.begin(); it != medoids.end(); it++){
		vector<int> &medoid = it->second;
		int NstatesInMedoid = medoid.size();
		for(int j=1;j<NstatesInMedoid;j++){ // Start at 1 because first index in medoid is center
			localStateNodes[medoid[j]].minDiv = bignum;
			for(int k=0;k<Nclu;k++){
				double div = wJSdiv(localStateNodes[medoid[j]].stateId,localStateNodes[k].stateId);
				if(div < localStateNodes[medoid[j]].minDiv){
					sumMinDiv -= localStateNodes[medoid[j]].minDiv;
					localStateNodes[medoid[j]].minDiv = div;
					sumMinDiv += localStateNodes[medoid[j]].minDiv;
					localStateNodes[medoid[j]].minCenterStateId = localStateNodes[k].stateId;
				}
			}	
		}
	}

	return sumMinDiv;

}

double StateNetwork::findCenters(vector<LocalStateNode> &localStateNodes){
	// Modifies the order of localStateNodes such thar the fist Nclu will be the centers.
	// Also, all elements will contain the stateId it is closest to.

	int NPstateNodes = localStateNodes.size();
	int Ncenters = 0;
	double sumMinDiv = bignum*NPstateNodes; // Because minDiv is set to bignum for all state nodes

	// Find random state node in physical node as first center
	std::uniform_int_distribution<int> randInt(0,NPstateNodes-1);
	int firstCenterIndex = randInt(mtRand);
	localStateNodes[firstCenterIndex].minCenterStateId = localStateNodes[firstCenterIndex].stateId;
	sumMinDiv -= localStateNodes[firstCenterIndex].minDiv;
	// Put the center in first non-center position (Ncenters = 0) by swapping elements
	swap(localStateNodes[Ncenters],localStateNodes[firstCenterIndex]);
	Ncenters++;


	// Find Nclu-1 more centers based on the k++ algorithm
	while(Ncenters < Nclu){
		int lastClusterId = localStateNodes[Ncenters-1].stateId;
		for(int i=Ncenters;i<NPstateNodes;i++){
			double div = wJSdiv(localStateNodes[i].stateId,lastClusterId);
			if(div < localStateNodes[i].minDiv){
				// Found new minimum divergence to center
				sumMinDiv -= localStateNodes[i].minDiv;
				sumMinDiv += div;
				localStateNodes[i].minDiv = div;
				localStateNodes[i].minCenterStateId = lastClusterId;
			}
		}
		// Pick new center proportional to minimum divergence
		std::uniform_int_distribution<double> randDouble(0.0,sumMinDiv);
		double randMinDivSum = randDouble(mtRand);
		double minDivSum = 0.0;
		int newCenterIndex = Ncenters;
		for(int i=Ncenters;i<NPstateNodes;i++){
			minDivSum += localStateNodes[i].minDiv;
			if(minDivSum > randMinDivSum){
				newCenterIndex = i;
				break;
			}
		}
		localStateNodes[newCenterIndex].minCenterStateId = localStateNodes[newCenterIndex].stateId;
		sumMinDiv -= localStateNodes[newCenterIndex].minDiv;
		// Put the center in first non-center position by swapping elements
		swap(localStateNodes[Ncenters],localStateNodes[newCenterIndex]);
		Ncenters++;
	}

	// Check if last center gives minimum divergence for some state nodes
	int lastClusterId = localStateNodes[Ncenters-1].stateId;
		for(int i=Ncenters;i<NPstateNodes;i++){
			double div = wJSdiv(localStateNodes[i].stateId,lastClusterId);
			if(div < localStateNodes[i].minDiv){
				// Found new minimum divergence to center
				sumMinDiv -= localStateNodes[i].minDiv;
				sumMinDiv += div;
				localStateNodes[i].minDiv = div;
				localStateNodes[i].minCenterStateId = lastClusterId;
			}
		}

		return sumMinDiv;

}

void StateNetwork::performLumping(vector<LocalStateNode> &localStateNodes){

	int NPstateNodes = localStateNodes.size();
	// Update stateNodes to reflect the lumping

	// // Validation
	// unordered_set<int> c;
	// for(int i=0;i<Nclu;i++)
	// 	c.insert(localStateNodes[i].stateId);
	// for(int i=Nclu;i<NPstateNodes;i++){
	// 	unordered_set<int>::iterator it = c.find(localStateNodes[i].minCenterStateId);
	// 	if(it == c.end())
	// 		cout << ":::::::::+++++++ ERROR for pos " << i << " " << localStateNodes[i].minDiv << " " << localStateNodes[i].minCenterStateId << endl;
	// }

	for(int i=Nclu;i<NPstateNodes;i++){


		// Only lump non-centers; first Nclu elements contain centers.

		StateNode &lumpedStateNode = stateNodes[localStateNodes[i].minCenterStateId];
		StateNode &lumpingStateNode = stateNodes[localStateNodes[i].stateId];
		// Add context to lumped state node
		lumpedStateNode.contexts.insert(lumpedStateNode.contexts.begin(),lumpingStateNode.contexts.begin(),lumpingStateNode.contexts.end());
		// Add links to lumped state node
		Nlinks  -= lumpedStateNode.links.size() + lumpingStateNode.links.size(); // To update the global number of links
		for(map<int,double>::iterator link_it = lumpingStateNode.links.begin(); link_it != lumpingStateNode.links.end(); link_it++){
			lumpedStateNode.links[link_it->first] += link_it->second;
		}
		Nlinks += lumpedStateNode.links.size(); // To update the global number of links

		lumpedStateNode.outWeight += lumpingStateNode.outWeight;

		// Update state id of lumping state node to point to lumped state node and make it inactive
		lumpingStateNode.updatedStateId = lumpedStateNode.stateId;
		lumpingStateNode.active = false;
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
			// The first Nclu elements will be centers
			vector<LocalStateNode> localStateNodes = vector<LocalStateNode>(NPstateNodes);
			for(int i=0;i<NPstateNodes;i++){
				localStateNodes[i].stateId = physNode.stateNodeIndices[i];
			}

			// Initialize vector to store centers and stateId of closest center
			double oldSumMinDiff = 0.0;
			double sumMinDiv = findCenters(localStateNodes);

			// Update centers as long as total distance to median changes more than threshold
			int Nupdates = 0;
			while( (fabs(sumMinDiv-oldSumMinDiff) > threshold) && (Nupdates < 10) ){
				swap(oldSumMinDiff,sumMinDiv);
				sumMinDiv = updateCenters(localStateNodes);
				Nupdates++;
			}

			// Free cached divergences
			cachedWJSdiv = unordered_map<pair<int,int>,double,pairhash>();

			// Perform the lumping and update stateNodes
			performLumping(localStateNodes);
			Nlumpings += NPstateNodes - Nclu;
			NstateNodes -= NPstateNodes - Nclu;

		}
		else if(NPstateNodes == 0){
			NphysDanglings++;
		}

		Nprocessed++;
		cout << "\r-->Lumped " << Nlumpings << " state nodes in " << Nprocessed << "/" << NphysNodes << " physical nodes.               ";
	}
	cout << endl << "-->Updating state node ids" << endl;

	// Update stateIds
	// First all active state nodes that other state nodes have lumped to
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(stateNode.active){
			stateNodeIdMapping[stateNode.stateId] = updatedStateId;
			stateNode.updatedStateId = updatedStateId;
			updatedStateId++;
		}
	}
	// Then all inactive state nodes that have lumped to other state nodes
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(!stateNode.active){
			stateNodeIdMapping[stateNode.stateId] = stateNodeIdMapping[stateNode.updatedStateId];
		}
	}

}

bool StateNetwork::readLines(string &line,vector<string> &lines){
	
	while(getline(ifs,line)){
		if(line[0] == '*'){
			return true;
		}
		else if(line[0] != '=' && line[0] != '#'){
			lines.push_back(line);
		}
	}

	return false; // Reached end of file
}

bool StateNetwork::loadStateNetworkBatch(){

	vector<string> stateLines;
	vector<string> linkLines;
	vector<string> contextLines;
	bool readStates = false;
	bool readLinks = false;
	bool readContexts = false;
	string buf;
	istringstream ss;

	// ************************* Read statenetwork batch ************************* //
	
	// Read until next data label. Return false if no more data labels
	if(keepReading){
		cout << "Reading statenetwork, batch " << Nbatches+1 << ":" << endl;
		if(line[0] != '*'){
			while(getline(ifs,line)){
				if(line[0] == '*')
					break;
				size_t foundTotSize = line.find("# Total weight: ");
				if(foundTotSize != string::npos)
					totWeight = atof(line.substr(foundTotSize+16).c_str());
			}
		}
	}
	else{
		cout << "-->No more statenetwork batches to read." << endl;
		return false;
	}

	while(!readStates || !readLinks || !readContexts){

		ss.clear();
		ss.str(line);
		ss >> buf;
		if(!readStates && buf == "*States"){
			cout << "-->Reading states..." << flush;
			readStates = true;
			keepReading = readLines(line,stateLines);
			NstateNodes = stateLines.size();
			cout << "found " << NstateNodes << " states." << endl;
		}
		else if(!readLinks && buf == "*Links"){
			cout << "-->Reading links..." << flush;
			readLinks = true;
			keepReading = readLines(line,linkLines);
			Nlinks = linkLines.size();
			cout << "found " << Nlinks << " links." << endl;
		}
		else if(!readContexts && buf == "*Contexts"){
			cout << "-->Reading contexts..." << flush;
			readContexts = true;
			keepReading = readLines(line,contextLines);
			Ncontexts = contextLines.size();
			cout << "found " << Ncontexts << " contexts." << endl;
		}
		else{
			cout << "Expected *States, *Links, or *Contexts, but found " << buf << " exiting..." << endl;
			exit(-1);
		}
	}

	// ************************* Process statenetwork batch ************************* //
	Nbatches++;
	cout << "Processing statenetwork, batch " << Nbatches << ":" << endl;

	//Process states
	cout << "-->Processing " << NstateNodes  << " state nodes..." << flush;
	for(int i=0;i<NstateNodes;i++){
		ss.clear();
		ss.str(stateLines[i]);
		ss >> buf;
		int stateId = atoi(buf.c_str());
		ss >> buf;
		int physId = atoi(buf.c_str());
	  ss >> buf;
	  double outWeight = atof(buf.c_str());
	  weight += outWeight;
		if(outWeight > epsilon)
			physNodes[physId].stateNodeIndices.push_back(stateId);
		else{
			physNodes[physId].stateNodeDanglingIndices.push_back(stateId);
			Ndanglings++;
		}
		stateNodes[stateId] = StateNode(stateId,physId,outWeight);
	}
	NphysNodes = physNodes.size();
	cout << "found " << Ndanglings << " dangling state nodes in " << NphysNodes << " physical nodes, done!" << endl;

	// Process links 
	cout << "-->Processing " << Nlinks  << " links..." << flush;
	for(int i=0;i<Nlinks;i++){
		ss.clear();
		ss.str(linkLines[i]);
		ss >> buf;
		int source = atoi(buf.c_str());
		ss >> buf;
		int target = atoi(buf.c_str());
		ss >> buf;
		double linkWeight = atof(buf.c_str());
		stateNodes[source].links[target] += linkWeight;
	}
 	cout << "done!" << endl;

	// Process contexts
	cout << "-->Processing " << Ncontexts  << " contexts..." << flush;
	for(int i=0;i<Ncontexts;i++){
		ss.clear();
		ss.str(contextLines[i]);
		ss >> buf;
		int stateNodeId = atoi(buf.c_str());
		string context = contextLines[i].substr(buf.length()+1);
		stateNodes[stateNodeId].contexts.push_back(context);
	}
	cout << "done!" << endl;

	// // Validate out-weights
 // 	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
	// StateNode &stateNode = it->second;
	// 	double w = 0.0;
	// 	for(map<int,double>::iterator it_link = stateNode.links.begin(); it_link != stateNode.links.end(); it_link++){
	// 		w += it_link->second;
	// 	}
	// 	if((w < (stateNode.outWeight-epsilon)) || (w > (stateNode.outWeight+epsilon))){
	// 		cout << setprecision(15) << "::::::::::: Warning: out-weight does not match link weights for state node " << stateNode.stateId << ": " << stateNode.outWeight << " vs " << w << " " << stateNode.links.size() << ", updating. :::::::::::" << endl;
	// 		stateNode.outWeight = w;
	// 	}
	// }

 	return true;

}


void StateNetwork::printStateNetworkBatch(){

	cout << "Writing temporary results:" << endl;

  my_ofstream ofs;
  if(batchOutput){
  	if(Nbatches == 1){ // Start with empty file for first batch
			ofs.open(tmpOutFileName.c_str());
		}
		else{ // Append to existing file
			ofs.open(tmpOutFileName.c_str(),ofstream::app);
		}
		ofs << "===== " << Nbatches << " =====\n";
		ofs << "*States\n";
		ofs << "#stateId ==> (physicalId, outWeight)\n";
  }
  else{
  	 if(Nbatches == 1){ // Start with empty file for first batch
			ofs.open(tmpOutFileNameStates.c_str());
		}
		else{ // Append to existing file
			ofs.open(tmpOutFileNameStates.c_str(),ofstream::app);
		}
  }

	cout << "-->Writing " << NstateNodes << " state nodes..." << flush;
	// To order state nodes by id
	map<int,int> orderedStateNodeIds;
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(stateNode.active)
			orderedStateNodeIds[stateNode.updatedStateId] = stateNode.stateId;
	}

	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){
		StateNode &stateNode = stateNodes[it->second];
		ofs << stateNode.stateId << " " << stateNode.physId << " " << stateNode.outWeight << "\n";
	}
	cout << "done!" << endl;

	if(batchOutput){
		ofs << "*Links\n";
		ofs << "#(source target) ==> weight\n";
	}
	else{
		ofs.close();
		if(Nbatches == 1){ // Start with empty file for first batch
			ofs.open(tmpOutFileNameLinks.c_str());
		}
		else{ // Append to existing file
			ofs.open(tmpOutFileNameLinks.c_str(),ofstream::app);
		}
	}
	cout << "-->Writing " << Nlinks << " links..." << flush;
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){	
		StateNode &stateNode = stateNodes[it->second];
		for(map<int,double>::iterator it_link = stateNode.links.begin(); it_link != stateNode.links.end(); it_link++){
				ofs << stateNode.stateId << " " << it_link->first << " " << it_link->second << "\n";
		}
	}
	cout << "done!" << endl;

	if(batchOutput){
		ofs << "*Contexts \n";
		ofs << "#stateId <== (physicalId priorId [history...])\n";
	}
	else{
		ofs.close();
		if(Nbatches == 1){ // Start with empty file for first batch
			ofs.open(tmpOutFileNameContexts.c_str());
		}
		else{ // Append to existing file
			ofs.open(tmpOutFileNameContexts.c_str(),ofstream::app);
		}
	}
	cout << "-->Writing " << Ncontexts << " contexts..." << flush;
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){	
		StateNode &stateNode = stateNodes[it->second];
		for(vector<string>::iterator it_context = stateNode.contexts.begin(); it_context != stateNode.contexts.end(); it_context++){
			ofs << stateNode.stateId << " " << (*it_context) << "\n";
		}
	}
	cout << "done!" << endl;

}

void StateNetwork::printStateNetwork(){

	entropyRate = calcEntropyRate();

  my_ofstream ofs;
  ofs.open(outFileName.c_str());

	cout << "No more batches, writing results to " << outFileName << ":" << endl;
	cout << "-->Writing header comments..." << flush;
  ofs << "# Number of physical nodes: " << NphysNodes << "\n";
  ofs << "# Number of state nodes: " << NstateNodes << "\n";
  ofs << "# Number of dangling physical (and state) nodes: " << NphysDanglings << "\n";
  ofs << "# Number of links: " << Nlinks << "\n";
  ofs << "# Number of contexts: " << Ncontexts << "\n";
  ofs << "# Total weight: " << weight << "\n";
  ofs << "# Entropy rate: " << entropyRate << "\n";
	cout << "done!" << endl;

	cout << "-->Writing " << NstateNodes << " state nodes..." << flush;
	// To order state nodes by id
	map<int,int> orderedStateNodeIds;
	for(unordered_map<int,StateNode>::iterator it = stateNodes.begin(); it != stateNodes.end(); it++){
		StateNode &stateNode = it->second;
		if(stateNode.active)
			orderedStateNodeIds[stateNode.updatedStateId] = stateNode.stateId;
	}

	ofs << "*States\n";
	ofs << "#stateId ==> (physicalId, outWeight)\n";
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){
		StateNode &stateNode = stateNodes[it->second];
		ofs << stateNodeIdMapping[stateNode.stateId] << " " << stateNode.physId << " " << stateNode.outWeight << "\n";
	}
	cout << "done!" << endl;

	cout << "-->Writing " << Nlinks << " links..." << flush;
	ofs << "*Links\n";
	ofs << "#(source target) ==> weight\n";
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){	
		StateNode &stateNode = stateNodes[it->second];
		int source = stateNodeIdMapping[stateNode.stateId];
		// Remove link redundance from lumped targets
		unordered_map<int,double> lumpedLinks;
		for(map<int,double>::iterator it_link = stateNode.links.begin(); it_link != stateNode.links.end(); it_link++)
			lumpedLinks[stateNodeIdMapping[it_link->first]] += it_link->second;
		for(unordered_map<int,double>::iterator it_link = lumpedLinks.begin(); it_link != lumpedLinks.end(); it_link++){
			ofs << source << " " << it_link->first << " " << it_link->second << "\n";
		}
	}
	cout << "done!" << endl;

	cout << "-->Writing " << Ncontexts << " contexts..." << flush;
	ofs << "*Contexts \n";
	ofs << "#stateId <== (physicalId priorId [history...])\n";
	for(map<int,int>::iterator it = orderedStateNodeIds.begin(); it != orderedStateNodeIds.end(); it++){	
		StateNode &stateNode = stateNodes[it->second];
		for(vector<string>::iterator it_context = stateNode.contexts.begin(); it_context != stateNode.contexts.end(); it_context++){
			ofs << stateNodeIdMapping[stateNode.stateId] << " " << (*it_context) << "\n";
		}
	}
	cout << "done!" << endl;

}

void StateNetwork::concludeBatch(){

	cout << "Concluding batch:" << endl;

	entropyRate += calcEntropyRate();
	accumWeight += weight;
	totNphysNodes += NphysNodes;
	totNstateNodes += NstateNodes;
	totNlinks += Nlinks;
	totNdanglings += Ndanglings;
	totNcontexts += Ncontexts;
	totNphysDanglings += NphysDanglings;
	weight = 0.0;
	NphysNodes = 0;
	NstateNodes = 0;
	Nlinks = 0;
	Ndanglings = 0;
	Ncontexts = 0;
	NphysDanglings = 0;

	cout << "-->Current estimate of the entropy rate: " << entropyRate*totWeight/accumWeight << endl;

	completeStateNodeIdMapping.insert(stateNodeIdMapping.begin(),stateNodeIdMapping.end());
	stateNodeIdMapping.clear();
	physNodes.clear();
	stateNodes.clear();
	cachedWJSdiv.clear();

}

void StateNetwork::compileBatches(){


  my_ofstream ofs;
  ofs.open(outFileName);
  string buf;
	istringstream ss;
	bool writeStates = false;
	bool writeLinks = false;
	bool writeContexts = false;
	int batchNr = 1;

	cout << "Writing final results to " << outFileName << ":" << endl;
  
  cout << "-->Writing header comments..." << flush;
  ofs << "# Number of physical nodes: " << totNphysNodes << "\n";
  ofs << "# Number of state nodes: " << totNstateNodes << "\n";
  ofs << "# Number of dangling physical (and state) nodes: " << totNphysDanglings << "\n";
  ofs << "# Number of links: " << totNlinks << "\n";
  ofs << "# Number of contexts: " << totNcontexts << "\n";
  ofs << "# Total weight: " << totWeight << "\n";
  ofs << "# Entropy rate: " << entropyRate << "\n";
	cout << "done!" << endl;

	cout << "-->Relabeling and writing " << totNstateNodes << " state nodes, " << totNlinks << " links, and " << totNcontexts << " contexts:" << endl;

	if(batchOutput){

		ifstream ifs_tmp(tmpOutFileName.c_str());
	
		// Copy lines directly until data format
		while(getline(ifs_tmp,line)){
			if(line[0] == '*'){
				break;	
			}
			else if(line[0] == '='){
				ofs << "=== " << batchNr << "/" << Nbatches << " ===\n";
			}
			else{
				ofs << line << "\n";
			}
		}
		while(!ifs_tmp.eof()){
	
			if(!writeStates && !writeLinks && !writeContexts){
				cout << "-->Batch " << batchNr << "/" << Nbatches << endl;
			}
			ofs << line << "\n";
			ss.clear();
			ss.str(line);
			ss >> buf;
			if(buf == "*States"){
				cout << "-->Writing state nodes..." << flush;
				writeStates = true;
				WriteMode writeMode = STATENODES;
				writeLines(ifs_tmp,ofs,writeMode,line,batchNr);
			}
			else if(buf == "*Links"){
				cout << "-->Writing links..." << flush;
				writeLinks = true;
				WriteMode writeMode = LINKS;
				writeLines(ifs_tmp,ofs,writeMode,line,batchNr);
			}
			else if(buf == "*Contexts"){
				cout << "-->Writing contexts..." << flush;
				writeContexts = true;
				WriteMode writeMode = CONTEXTS;
				writeLines(ifs_tmp,ofs,writeMode,line,batchNr);
			}
			else{
				cout << "Failed on line: " << line << endl;
			}
			cout << "done!" << endl;
			if(writeStates && writeLinks && writeContexts){
				writeStates = false;
				writeLinks = false;
				writeContexts = false;
				batchNr++;
			}
		}
		remove( tmpOutFileName.c_str() );
	}
	else{

		cout << "-->Writing " << totNstateNodes << " state nodes..." << flush;
		ofs << "*States\n";
		ofs << "#stateId ==> (physicalId, outWeight)\n";
		ifstream ifs_tmp(tmpOutFileNameStates.c_str());
		WriteMode writeMode = STATENODES;
		writeLines(ifs_tmp,ofs,writeMode,line);	
		ifs_tmp.close();
		cout << "done!" << endl;


		cout << "-->Writing " << totNlinks << " links..." << flush;
		ofs << "*Links\n";
		ofs << "#(source target) ==> weight\n";
		ifs_tmp.open(tmpOutFileNameLinks.c_str());
		writeMode = LINKS;
		writeLines(ifs_tmp,ofs,writeMode,line);	
		ifs_tmp.close();
		cout << "done!" << endl;

		cout << "-->Writing " << totNcontexts << " contexts..." << flush;
		ofs << "*Contexts \n";
		ofs << "#stateId <== (physicalId priorId [history...])\n";		
		ifs_tmp.open(tmpOutFileNameContexts.c_str());		
		writeMode = CONTEXTS;
		writeLines(ifs_tmp,ofs,writeMode,line);	
		ifs_tmp.close();
		cout << "done!" << endl;

		remove( tmpOutFileNameStates.c_str() );
		remove( tmpOutFileNameLinks.c_str() );
		remove( tmpOutFileNameContexts.c_str() );
	}

}

void StateNetwork::writeLines(ifstream &ifs_tmp, ofstream &ofs, WriteMode &writeMode, string &line,int &batchNr){

	string buf;
	istringstream ss;

	unordered_map<pair<int,int>,double,pairhash> aggregatedLinks;

	while(getline(ifs_tmp,line)){
		if(line[0] != '*'){
			if(line[0] != '=' && line[0] != '#'){
				ss.clear();
				ss.str(line);
				ss >> buf;
				if(writeMode == STATENODES){
					int stateId = atoi(buf.c_str());
					ss >> buf;
					int physId = atoi(buf.c_str());
	 				ss >> buf;
	 				double outWeight = atof(buf.c_str());
					ofs << completeStateNodeIdMapping[stateId] << " " << physId << " " << outWeight << "\n";
				}
				else if(writeMode == LINKS){
					int source = atoi(buf.c_str());
					ss >> buf;
					int target = atoi(buf.c_str());
					ss >> buf;
					double linkWeight = atof(buf.c_str());
					ofs << completeStateNodeIdMapping[source] << " " << completeStateNodeIdMapping[target] << " " << linkWeight << "\n";
					aggregatedLinks[make_pair(source,target)] += linkWeight;				
				}
				else if(writeMode == CONTEXTS){
					int stateNodeId = atoi(buf.c_str());
					string context = line.substr(buf.length()+1);
					ofs << completeStateNodeIdMapping[stateNodeId] << " " << context << "\n";
				}
			}
			else if(line[0] == '='){
				ofs << "=== " << batchNr+1 << "/" << Nbatches << " ===\n";
			}
			else{
				ofs << line << "\n";
			}
		}
		else{
			break;
		}
	}

	if(writeMode == LINKS ){
		for(unordered_map<pair<int,int>,double,pairhash>::iterator link_it = aggregatedLinks.begin(); link_it != aggregatedLinks.end(); link_it++){
			ofs << link_it->first.first << " " << link_it->first.second << " " << link_it->second << "\n";
		}
	}
}

void StateNetwork::writeLines(ifstream &ifs_tmp, ofstream &ofs, WriteMode &writeMode, string &line){

	string buf;
	istringstream ss;

	while(getline(ifs_tmp,line)){
		ss.clear();
		ss.str(line);
		ss >> buf;
		if(writeMode == STATENODES){
			int stateId = atoi(buf.c_str());
			ss >> buf;
			int physId = atoi(buf.c_str());
	 		ss >> buf;
	 		double outWeight = atof(buf.c_str());
			ofs << completeStateNodeIdMapping[stateId] << " " << physId << " " << outWeight << "\n";
		}
		else if(writeMode == LINKS){
			int source = atoi(buf.c_str());
			ss >> buf;
			int target = atoi(buf.c_str());
			ss >> buf;
			double linkWeight = atof(buf.c_str());
			ofs << completeStateNodeIdMapping[source] << " " << completeStateNodeIdMapping[target] << " " << linkWeight << "\n";					
		}
		else if(writeMode == CONTEXTS){
			int stateNodeId = atoi(buf.c_str());
			string context = line.substr(buf.length()+1);
			ofs << completeStateNodeIdMapping[stateNodeId] << " " << context << "\n";
		}
	}
}
