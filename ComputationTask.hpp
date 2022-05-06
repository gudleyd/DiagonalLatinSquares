//
// Created by Ivan Lebedev on 4/25/22.
//

#pragma once

#include "DLS.hpp"

#include <vector>
#include <cstring>
#include <fstream>


class ComputationTask {
public:
    ComputationTask(const char* inputFileName, const char* outputFileName):
        taskData_(ReadTask(inputFileName)),
        initialState_(GenerateInitialState(taskData_)),
        allSquares_(DiagonalLatinSquareGenerator::GenerateAllPossibleSquares(initialState_)) {

        out_.open(outputFileName, std::fstream::out | std::fstream::app | std::fstream::binary);
        outputFileName_ = outputFileName;
    }

    ~ComputationTask() {
        out_.flush();
        out_.close();
    }

    void GoToCheckpoint(const char* checkpointFileName) {
        std::fstream checkpoint;
        checkpoint.open(checkpointFileName, std::ios::in | std::fstream::binary);
        if (checkpoint.fail() || checkpoint.tellp() == 0) return;
        uint64_t nchars;
        checkpoint >> taskData_.processingPosition >> nchars;
        // open output file and truncate it to nchars
        out_ = std::fstream(outputFileName_, std::fstream::out | std::fstream::app | std::fstream::binary | std::fstream::trunc);
        out_.seekp(nchars, std::fstream::beg);
    }

    bool DoIteration() {
        if (taskData_.processingPosition >= allSquares_.size()) {
            return false;
        }
        auto currentSquare = allSquares_[taskData_.processingPosition];
        auto sublatins = DLS_SubrectangleFinder::CountSubLatin(currentSquare);
        auto loops = DLS_LoopFinder::CountLoops(currentSquare);
        out_ << currentSquare.ToString() << " " << sublatins << " " << loops.first <<  " " << loops.second << std::endl;
        return ++taskData_.processingPosition != allSquares_.size();
    }

    int MakeCheckpoint(const char* checkpointFileName) {
        out_.flush();
        std::fstream checkpoint;
        checkpoint.open(checkpointFileName, std::ios::out | std::fstream::binary);
        if (checkpoint.fail()) return -1;
        checkpoint << taskData_.processingPosition << " " << out_.tellp() << std::endl;
        checkpoint.close();
        return 0;
    }
    
    void Finish() {
        out_.flush();
    }

    double GetFractionDone() {
        return static_cast<double>(taskData_.processingPosition) / allSquares_.size();
    }

    struct TaskData {
        uint64_t size = 1;
        uint64_t positions[7] = {0};

        // checkpoint
        uint64_t processingPosition = 0;
    };

    static bool PossibleToInitialize(TaskData taskData) {
        try {
            GenerateInitialState(taskData);
        } catch (...) {
            return false;
        }
        return true;
    }

private:

    TaskData taskData_;
    DiagonalLatinSquareGenerationState initialState_;
    std::vector<DiagonalLatinSquareGenerationState> allSquares_;
    std::fstream out_;
    const char* outputFileName_;

private:

    static TaskData ReadTask(const char* inputFileName) {
        std::fstream in(inputFileName, std::ios::in);
        uint64_t size, pos[7];
        if (!(in >> size >> pos[0] >> pos[1] >> pos[2] >> pos[3] >> pos[4] >> pos[5] >> pos[6])) {
            Application::fatalError("Wrong input format");
        }
        TaskData taskData;
        taskData.size = size;
        std::memcpy(taskData.positions, pos, sizeof(pos));
        return taskData;
    }

    static DiagonalLatinSquareGenerationState GenerateInitialState(const TaskData& taskData) {
        DiagonalLatinSquareGenerationState state(taskData.size);
        auto mainDiagonal = SequenceGenerator::GenerateSequence<uint8_t>(taskData.size, taskData.positions[0]);
        for (uint8_t i = 0; i < state.Size(); ++i) {
            state.SetValue(i, i, mainDiagonal[i]);
        }

        std::vector<std::vector<uint8_t>> possibleValues;
        for (uint8_t i = 0; i < state.Size(); ++i) {
            possibleValues.emplace_back(state.GetPossibleValues(i, state.Size() - 1 - i));
        }

        auto antiDiagonal = SequenceGenerator::GenerateSequence(possibleValues, taskData.positions[1]);
        for (uint8_t i = 0; i < state.Size(); ++i) {
            state.SetValue(i, state.Size() - i - 1, antiDiagonal[i]);
        }

        GenerateRow(state, 0, taskData.positions[2]);
        GenerateColumn(state, 0, taskData.positions[3]);
        GenerateRow(state, 1, taskData.positions[4]);
        GenerateColumn(state, 1, taskData.positions[5]);
        GenerateRow(state, 2, taskData.positions[6]);

        return state;
    }

    static void GenerateRow(DiagonalLatinSquareGenerationState& state, uint8_t rowIndex, uint8_t sequenceIndex) {
        std::vector<std::vector<uint8_t>> possibleValues;
        for (uint8_t i = 0; i < state.Size(); ++i) {
            possibleValues.emplace_back(state.GetPossibleValues(rowIndex, i));
        }
//        std::cout << (int)rowIndex << " " << SequenceGenerator::CountSequences(possibleValues) << std::endl;
        auto column = SequenceGenerator::GenerateSequence(possibleValues, sequenceIndex);
        for (uint8_t i = 0; i < state.Size(); ++i) {
            state.SetValue(rowIndex, i, column[i]);
        }
    }

    static void GenerateColumn(DiagonalLatinSquareGenerationState& state, uint8_t columnIndex, uint8_t sequenceIndex) {
        std::vector<std::vector<uint8_t>> possibleValues;
        for (uint8_t i = 0; i < state.Size(); ++i) {
            possibleValues.emplace_back(state.GetPossibleValues(i, columnIndex));
        }
        auto column = SequenceGenerator::GenerateSequence(possibleValues, sequenceIndex);
        for (uint8_t i = 0; i < state.Size(); ++i) {
            state.SetValue(i, columnIndex, column[i]);
        }
    }
};