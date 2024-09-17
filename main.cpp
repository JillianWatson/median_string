#include <iostream> 
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <climits>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

using namespace std;

//define alphabet 
const vector<char> NT = {'A','C','G','T'};


// calc hamming distance between 2 strings
int hammingDistance(const string& str1, const string& str2) {
    if(str1.length() != str2.length()) {
        throw invalid_argument("k-mers must be equal length");
    }
    int dist = 0;
    for (size_t i =0; i < str1.length(); ++i) {
        if (str1[i] != str2[i]) {
                dist++;
        }
    }
     
    return dist;
}

// find minimum distance between 2 strings
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

// calculate total distance between 2 strings
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

// check that file contents are loaded correctly
void checkSequences(const vector<string>& sequences, int NT = 10){
    cout << "Checking first " << NT << " nucleotides of each sequence: " << endl;
    for (size_t i=0; i<sequences.size(); ++i) {
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

    // reached leaf node
    if (iter == K) {
        //cout << "evaluating kmer: " << currentStr << " with distance: " << bestDistance << endl;
        int distance = distanceTotal(currentStr, sequences);

        if (distance < bestDistance) {
            bestDistance = distance;
            bestStr = currentStr;
            //cout << "New best k-mer found: " << bestStr << " with distance: " << bestDistance << endl;
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
        cerr << "Please provide one multi-sequence input file." << endl;
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
    int K;
    cout << "Provide desired length of k-mer: ";
    cin >> K;

    if (K <= 0 || K > 10) {
        cerr << "Error: please provide k-mer length between 1 and 9." << endl;
        return 1;
    }

    checkSequences(sequences);
    for (size_t i=0; i< sequences.size(); ++i) {
        cout << "Sequence: " << i+1 << " length: " << sequences[i].length() << endl; 
    }

    string bestStr(K,'A');
    string currentStr(K,'A');
    int bestDistance = INT_MAX;

    cout << "Initial string: " << bestStr << "distance: " <<bestDistance << endl;
    
    cout << "Starting branch and bound algo with K = " << K << endl;
    branch_and_bound(sequences, currentStr, bestStr, bestDistance, 0, K);

    cout << "best string: " << bestStr << endl;
    cout << "distance: " << bestDistance << endl;

    return 0;
 }

