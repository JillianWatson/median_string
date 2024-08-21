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

//alphabet for creating k-mers
const vector<char> NT = {'A','C','G','T'};

// calculate hamming distance between 2 strings
int hammingDistance(const string& str1, const string& str2) {
    if(str1.length() != str2.length()) {
        throw invalid_argument("k-mers must be equal length");
    }
    int dist = 0;
    for (size_t i =0; i < str1.length(); i++) {
        if (str1[i] != str2[i]) {
                dist ++;
        }
    }
    return dist;
}
// calculate distance between kmer and a single sequence
int distanceToSequence(const string& kmer, const string& seq) {
    int minDist = numeric_limits<int>::max();
    for (size_t i=0; i <= seq.length() - kmer.length(); i++) {
        minDist = min(minDist, hammingDistance(kmer, seq.substr(i, kmer.length())));
    }
    return minDist;

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
    
    // grab user input; determine length of k-mer
    int K;
    cout << "Provide desired length of k-mer: ";
    cin >> K;

    if (K <= 0 || K > 10) {
        cerr << "Error: please provide k-mer length between 1 and 9." << endl;
        return 1;
    }

    for (size_t i=0; i< sequences.size(); ++i) {
        cout << "Sequence " << i+1 << " length: " << sequences[i].length() << endl; 
    }

    return 0;
}
