#include <iostream> 
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <climits>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <random>
#include <iterator>

using namespace std;


//define alphabet 
const vector<char> NT = {'A','C','G','T'};


// calc hamming distance between 2 strings
int hammingDistance(const string& str1, const string& str2) {
    if(str1.length() != str2.length()) {
        throw invalid_argument("strings must be of equal length");
    }
    int dist = 0;
    for (size_t i=0; i < str1.length(); ++i) {
        if (str1[i] != str2[i]) {
                dist++;
        }
    }
     
    return dist;
}


// helper function for distance calculation
int distanceToSequence(const string& kmer, const string& seq) {
    int minDist = numeric_limits<int>::max();
    for (size_t i=0; i <= seq.length() - kmer.length(); i++) {
        int dist = hammingDistance(kmer, seq.substr(i, kmer.length()));
        if (dist < minDist) {
            minDist = dist;
        }
    }
  
    return minDist;
}


// calculate distance between 2 strings
int distanceTotal(const string& kmer, const vector<string>& sequences) {
    int total = 0;
    //cout << "calculating total distance for k-mer: " << kmer << endl;
    for (const auto& seq : sequences) {
        int dist = distanceToSequence(kmer, seq);
        total += dist;
        //cout << "distance for kmer " << kmer << " : " << dist << endl;
    }
    //cout << "total dist: " << total << endl;
    return total;
}


// find a more suitable starting string at the cost of reproducability 
// builds the starting k-mer from sampling of input sequences, calculating the most frequent nucleotide at each position
// risk of introducing bias 
string HeuristicKmer(const vector<string>& sequences, int K) {
    array<int, 4> countNT = {0,0,0,0};
    int counter = 0;

    // seeding random generator
    random_device rd;
    mt19937 gen(rd());
    int sampleSize = min(10000, static_cast<int>(sequences.size() * K));

    // sampling loop. select random sequence, seqIter, and random start position, posIter, for NT counting
    for (int i = 0; i < sampleSize; ++i) {
        int seqIter = gen() % sequences.size();
        int posIter = gen() % (sequences[seqIter].length() - K + 1);

        char c = sequences[seqIter][posIter];
        // if A,C,G,T found in sequence, increment counter for that NT. 
        auto j = find(NT.begin(), NT.end(), c);
        if (j != NT.end()) {
            countNT[distance(NT.begin(), j)]++;
            counter++;
        }
    }

    // build starting kmer
    string startKmer(K, 'A');
    for (int i=0; i < K; i++){
        int maxCount = 0;
        int maxIndex = 0;
        for (int j=0; j < NT.size(); j++){
            if (countNT[j] > maxCount) {
                maxCount = countNT[j];
                maxIndex = j;    
            }
        }
        startKmer[i] = NT[maxIndex];
        // zero NT counts after selection of each position
        countNT[maxIndex] = 0;
    }

    return startKmer;
            
}


// check that file contents are loaded correctly
void checkSequences(const vector<string>& sequences, int NT = 10){
    cout << "Checking first " << NT << " nucleotides of each sequence: " << endl;
    for (size_t i=0; i < sequences.size(); ++i) {
        cout << "Sequence " << i + 1 << ": ";
        if (sequences[i].length() < NT) {
            cout << sequences[i] << " (full sequence; total length: " << sequences[i].length() << ")";
        } else {
            cout << sequences[i].substr(0, NT) << " (partial sequence; total length: " << sequences[i].length() << ")";
        }
        cout << endl;    
    }
    cout << endl;

}


void branch_and_bound(const vector<string>& sequences, string& currentStr, string& bestStr, int& bestDistance, int iter, int K) {

    int currentDistance = distanceTotal(currentStr.substr(0, iter), sequences);

    if (currentDistance >= bestDistance) {
        return;
    }

    // reached leaf node
    if (iter == K) {
        if (currentDistance < bestDistance) {
            bestDistance = currentDistance;
            bestStr = currentStr;
        }
        return;
    }

    for (char nucleotide : NT) {
        currentStr[iter] = nucleotide;
        //cout << "Exploring current nucleotide " << nucleotide << " at position " << iter << endl;
        branch_and_bound(sequences, currentStr, bestStr, bestDistance, iter + 1, K);    
    }

}


int main(int argc, char* argv[]) {

    if (argc !=2) {
        cerr << "Please provide one multi-fasta input file." << endl;
        return 1;
    }

    ifstream inputFile(argv[1]);

    if (!inputFile) {
        cerr << "Error: unable to open input file." << endl;
        return 1;
    }

    string line, sequence;
    vector<string> sequences;

    while (getline(inputFile, line)) {
        // check for beginning of fasta sequence
        if (line[0] == '>') {
            // if sequence exists, add to sequences vector. Then clear sequence string to parse next seq string.
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                sequence.clear();
            }
        } else {
            sequence += line;
        }
    }
    if (!sequence.empty()) {
        sequences.push_back(sequence);
    }
    
    inputFile.close();
    
    // grab user input; determine length of desired k-mer
    // when K <= 4, premature pruning occurs. given small search space, brute force may be used.
    // when K > 10, performance decreases. Implement optimization techniques like multithreading, bit representation...
    int K;
    cout << "Provide desired length of k-mer: ";
    cin >> K;

    if (K <= 4 || K > 10) {
        cerr << "Error: please provide k-mer length between 4 and 9." << endl;
        return 1;
    }

    //checkSequences(sequences);
    for (size_t i=0; i< sequences.size(); ++i) {
        cout << "Sequence: " << i+1 << " length: " << sequences[i].length() << endl; 
    }
    
    cout << endl;

    // for naive branch and bound 
    string bestStr(K,'A');
    string currentStr = bestStr;
    int bestDistance = INT_MAX;

    // for heuristic b&b
    string heurBestStr = HeuristicKmer(sequences, K);
    string heurCurrentStr = heurBestStr;
    int heurBestDistance = distanceTotal(heurCurrentStr, sequences);

    cout << "Heuristic initial string: " << heurBestStr << " with start distance: " << heurBestDistance << endl;
    cout << "Starting heuristic branch and bound with K = " << K << endl;
    branch_and_bound(sequences, heurCurrentStr, heurBestStr, heurBestDistance, 0, K);

    cout << endl; 
    cout << "heuristic final best string: " << heurBestStr << " with final distance: " << heurBestDistance << endl;
    cout << endl;
    cout << endl;

    cout << "Naive initial string: " << bestStr << " with start distance: " << bestDistance << endl;
    cout << "Starting naive branch and bound algo with K = " << K << endl;
    branch_and_bound(sequences, currentStr, bestStr, bestDistance, 0, K);
    
    cout << endl;
    cout << "naive final best string: " << bestStr << " with final distance: " << bestDistance << endl;
    cout << endl;
    return 0;
 }

