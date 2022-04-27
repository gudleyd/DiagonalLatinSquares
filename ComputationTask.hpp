//
// Created by Ivan Lebedev on 4/25/22.
//

#pragma once

#include "DLS.hpp"

#include <fstream>


class ComputationTask {
public:
    static void Run(const std::string& inputFileName, const std::string& outputFileName) {
        auto taskData = ReadTask(inputFileName);
        auto initialState = GenerateInitialState(taskData);
        auto allSquares = DiagonalLatinSquareGenerator::GenerateAllPossibleSquares(initialState);
        // measure time for code to run
        auto start = std::chrono::high_resolution_clock::now();
        for (const auto &square: allSquares) {
            std::cout << DLS_SubrectangleFinder::CountSubLatin(square) << std::endl;
            std::cout << DLS_LoopFinder::CountLoops(square).first << std::endl;
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;
    }

private:

    struct TaskData {
        uint64_t size = 1;
        uint64_t positions[7] = {0};
    };

    static TaskData ReadTask(const std::string& inputFileName) {
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