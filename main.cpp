#include "DLS.hpp"

void PerformanceTest();

// print std::set operator ostream
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s)
{
    os << "{";
    for (const auto& x : s)
    {
        os << x << ", ";
    }
    os << "}";
    return os;
}

void kek() {
    int n = 10;
    std::vector<std::set<int>> set(5);
    for (int iii = 0; iii < 100; ++iii) {
        std::cout << iii << std::endl;
        DiagonalLatinSquareGenerationState state(n);
        auto mainDiagonal = SequenceGenerator::GenerateSequence<int>(10, rand() % 3628800);
        for (int i = 0; i < n; ++i) {
            state.SetValue(i, i, mainDiagonal[i]);
        }

        std::vector<std::vector<uint8_t>> possibleValues;
        for (int i = 0; i < n; ++i) {
            possibleValues.push_back(state.GetPossibleValues(i, n - 1 - i));
        }

        set[0].insert(SequenceGenerator::CountSequences(possibleValues));
        auto antiDiagonal = SequenceGenerator::GenerateSequence(possibleValues, 109);
        for (int i = 0; i < n; ++i) {
            state.SetValue(i, n - i - 1, antiDiagonal[i]);
        }

        possibleValues.clear();
        for (int i = 0; i < n; ++i) {
            possibleValues.push_back(state.GetPossibleValues(0, i));
        }

        set[1].insert(SequenceGenerator::CountSequences(possibleValues));
        auto firstRow = SequenceGenerator::GenerateSequence(possibleValues, 100);
        for (int i = 0; i < n; ++i) {
            state.SetValue(0, i, firstRow[i]);
        }

        possibleValues.clear();
        for (int i = 0; i < n; ++i) {
            possibleValues.push_back(state.GetPossibleValues(i, 0));
        }

        set[2].insert(SequenceGenerator::CountSequences(possibleValues));
        auto firstColumn = SequenceGenerator::GenerateSequence(possibleValues, 100);
        for (int i = 0; i < n; ++i) {
            state.SetValue(i, 0, firstColumn[i]);
        }

        possibleValues.clear();
        for (int i = 0; i < n; ++i) {
            possibleValues.push_back(state.GetPossibleValues(1, i));
        }

        set[3].insert(SequenceGenerator::CountSequences(possibleValues));
        auto secondRow = SequenceGenerator::GenerateSequence(possibleValues, 100);
        for (int i = 0; i < n; ++i) {
            state.SetValue(1, i, secondRow[i]);
        }

        possibleValues.clear();
        for (int i = 0; i < n; ++i) {
            possibleValues.push_back(state.GetPossibleValues(i, 1));
        }

        set[4].insert(SequenceGenerator::CountSequences(possibleValues));
        auto secondColumn = SequenceGenerator::GenerateSequence(possibleValues, 100);
        for (int i = 0; i < n; ++i) {
            state.SetValue(i, 1, secondColumn[i]);
        }

//        std::cout << state << std::endl;

//        state.PrintMatrix();
//
//        // measure time
//        auto start = std::chrono::high_resolution_clock::now();
//        auto allSquares = DiagonalLatinSquareGenerator::GenerateAllPossibleSquares(state);
//        std::cout << "All squares: " << allSquares.size() << std::endl;
//        auto end = std::chrono::high_resolution_clock::now();
//        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//        std::cout << "Time: " << duration.count() << " ms" << std::endl;
    }

    for (size_t i = 0; i < set.size(); ++i) {
        std::cout << "i: " << i << " count: " << set[i].size() << std::endl;
        std::cout << set[i] << std::endl;
    }
}

#include "ComputationTask.hpp"

int main() {
//    const std::string path = "kek_456";
//    std::ifstream result(path, std::ifstream::in);
//    std::string srName, seqno;
//    for (int j=path.size() - 1; j>=0 && path[j]!='_'; --j) seqno = path[j] + seqno;
//    uint64_t subrectangles, full_cycles, partial_cycles;
//    std::string currentTransaction;
//    int currentPos = 0;
//    while (result >> srName >> subrectangles >> full_cycles >> partial_cycles) {
//        currentTransaction += "INSERT OR IGNORE INTO RESULT(ID, SEQNO, POSINSEQ, SUBRECTANGLES, FULL_CYCLES, PARTIAL_CYCLES) VALUES ('" + srName + "', " + seqno + ", " + std::to_string(currentPos) + ", " + std::to_string(subrectangles) + ", " + std::to_string(full_cycles) + ", " + std::to_string(partial_cycles) + ");";
//        currentPos += 1;
//    }
//    std::cout << currentTransaction << std::endl;
//    return 0;
    auto start = std::chrono::high_resolution_clock::now();

    auto task = ComputationTask("in.txt", "out.txt");
    task.GoToCheckpoint("checkpoint.txt");
    while (task.DoIteration()) {
        std::cout << task.GetFractionDone() << std::endl;
        task.MakeCheckpoint("checkpoint.txt");
        task.GoToCheckpoint("checkpoint.txt");
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time: " << duration.count() << " ms" << std::endl;
    return 0;
//    kek();
    return 0;
//    std::cout << SequenceGenerator::CountSequences<int>(vv) << std::endl;
//    std::cout << SequenceGenerator::GenerateSequence<int>(11, 3000000) << std::endl;
    return 0;
    auto squares = DiagonalLatinSquareGenerator(7).GenerateSquares();
    auto startTime = std::chrono::high_resolution_clock::now();
    {
//        auto squares = DiagonalLatinSquareGenerator(8).GenerateSquares();
        uint32_t maximumFullCycles = 0, maximumPartialCycles = 0;
        DiagonalLatinSquare maximumFullCyclesSquare(squares.front().Size()), maximumPartialCyclesSquare(squares.front().Size());
        for (const auto& square : squares) {
            auto [fullCycles, partialCycles] = DLS_LoopFinder::CountLoops(square);
            if (maximumFullCycles < fullCycles) {
                maximumFullCycles = fullCycles;
                maximumFullCyclesSquare = square;
            }
            if (maximumPartialCycles < partialCycles) {
                maximumPartialCycles = partialCycles;
                maximumPartialCyclesSquare = square;
            }
        }
        std::cout << "Maximum full cycles: " << maximumFullCycles << std::endl;
        std::cout << "Square with maximum full cycles:\n" << maximumFullCyclesSquare << std::endl;
        std::cout << "Maximum partial cycles: " << maximumPartialCycles << std::endl;
        std::cout << "Square with maximum partial cycles:\n" << maximumPartialCyclesSquare << std::endl;
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count() << " ms" << std::endl;
    return 0;
}

void PerformanceTest() {
    int numCount = 256, repeatCount = 1000000;
    std::vector<int> vector;
    std::set<int> set;
    std::vector<std::set<int>> vectorSet;
    for (int i = 0; i < numCount; ++i) {
        vector.emplace_back(i);
        set.emplace(i);
    }
    for (int i = 0; i < repeatCount; ++i) {
        vectorSet.emplace_back(set);
    }

    // Measure Time
    auto startTime = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < repeatCount; ++k) {
        auto lvector = vector;
        for (int i = 0; i < numCount; ++i) {
            lvector.erase(std::find(lvector.begin(), lvector.end(), i));
        }
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Vector Time : " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms" << std::endl;

    startTime = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < repeatCount; ++k) {
        auto lset = std::move(vectorSet[k]);
        for (int i = 0; i < numCount; ++i) {
            lset.erase(i);
        }
    }
    endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Set Time : " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms" << std::endl;
}