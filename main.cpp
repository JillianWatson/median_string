#include <iostream> 
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <climits>
#include <cmath>

using namespace std;

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

    string readline, sequence;
    vector<string> sequences;

    while (getline(inputFile, readline)) {
        if (readline[0] == '>') {
            // if sequence exists, add to the sequences vector. Then clear the sequence string for next string.
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                sequence.clear();
            }
        } else {
            sequence += readline;
        }
    }
    if (!sequence.empty()) {
        sequences.push_back(sequence);
    }

    for (size_t i=0; i< sequences.size(); ++i) {
        cout << "Sequence" << i+1 << " length: " << sequences[i].length() << endl; 
    }

    return 0;
}
