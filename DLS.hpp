#pragma once

#include <cinttypes>
#include <random>
#include <iostream>
#include <new>
#include <list>
#include <set>
#include <utility>
#include <numeric>
#include <queue>
#include <type_traits>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <chrono>

class Application {
public:
    static void fatalError(const char* message) {
        std::flush(std::cout);
        std::cerr << message << std::endl;
        std::flush(std::cerr);
        std::exit(1);
    }
};

// Fast random number generator
class RandomGenerator {
public:
    static RandomGenerator& Shared() {
        static RandomGenerator generator;
        return generator;
    }

    explicit RandomGenerator(uint32_t seed = time(nullptr)) {
        generator = std::mt19937(seed);
    }

    void SetSeed(uint32_t seed) {
        generator.seed(seed);
    }

    uint32_t Generate() {
        return generator();
    }

    uint32_t Generate(uint32_t max) {
        return generator() % max;
    }

    template<typename T>
    T GetRandomValue(const std::vector<T>& values) {
        return values[Generate() % values.size()];
    }

    template<typename T>
    T GetRandomValue(const std::set<T>& values) {
        auto it = values.begin();
        for (uint32_t i = Generate() % values.size(); i > 0; --i) {
            ++it;
        }
        return *it;
    }

    RandomGenerator(RandomGenerator const&) = delete;
    void operator=(RandomGenerator const&) = delete;

private:
    std::mt19937 generator;
};

class DiagonalLatinSquare {
public:
    static uint8_t Empty;

    explicit DiagonalLatinSquare(uint8_t size) {
        square_ = std::vector<uint8_t>(size * size, Empty);
        size_ = size;
    }

    [[nodiscard]] static DiagonalLatinSquare FromString(const std::string& str) {
        auto sizeRoot = static_cast<uint8_t>(std::sqrt(str.size()));
        if (sizeRoot * sizeRoot != str.size() || (sizeRoot > 1 && sizeRoot < 4)) {
            Application::fatalError("Invalid string size");
        }

        DiagonalLatinSquare square(sizeRoot);
        for (size_t i = 0; i < str.size(); ++i) {
            if (str[i] == '.') {
                square(i / square.Size(), i % square.Size()) = Empty;
            } else if (isdigit(str[i])) {
                square(i / square.Size(), i % square.Size()) = str[i] - '0';
            } else {
                square(i / square.Size(), i % square.Size()) = toupper(str[i]) - 'A' - '0';
            }

        }
        return square;
    }

    [[nodiscard]] std::string ToString() const {
        std::string str(size_ * size_, '.');
        for (uint8_t i = 0; i < size_ * size_; ++i) {
            if (square_[i] == Empty) {
                str[i] = '.';
            } else if (square_[i] < 10) {
                str[i] = square_[i] + '0';
            } else {
                str[i] = square_[i] - 10 + 'A';
            }
        }
        return str;
    }

    [[nodiscard]] uint8_t Size() const {
        return size_;
    }

    [[nodiscard]] bool AreSizesEqual() const {
        return size_ * size_ == square_.size();
    }

    [[nodiscard]] uint8_t& operator()(uint8_t i, uint8_t j) {
        return square_[i * size_ + j];
    }

    [[nodiscard]] const uint8_t& operator()(uint8_t i, uint8_t j) const {
        return square_[i * size_ + j];
    }

    [[nodiscard]] bool IsValid() const {
        for (uint8_t i = 0; i < size_; ++i) {
            for (uint8_t j = 0; j < size_; ++j) {
                if (square_[i * size_ + j] == Empty) {
                    return false;
                }
            }
        }

        // Check each row contains unique elements
        for (uint8_t i = 0; i < size_; ++i) {
            std::set<uint8_t> row;
            for (uint8_t j = 0; j < size_; ++j)
                row.insert(square_[i * size_ + j]);
            if (row.size() != size_)
                return false;
        }

        // Check each column contains unique elements
        for (uint8_t j = 0; j < size_; ++j) {
            std::set<uint8_t> column;
            for (uint8_t i = 0; i < size_; ++i)
                column.insert(square_[i * size_ + j]);
            if (column.size() != size_)
                return false;
        }

        // Check main diagonal contains unique elements
        std::set<uint8_t> diagonal;
        for (uint8_t i = 0; i < size_; ++i)
            diagonal.insert(square_[i * size_ + i]);
        if (diagonal.size() != size_)
            return false;

        // Check anti-diagonal contains unique elements
        diagonal.clear();
        for (uint8_t i = 0; i < size_; ++i)
            diagonal.insert(square_[i * size_ + size_ - i - 1]);
        if (diagonal.size() != size_)
            return false;

        return true;
    }

protected:
    std::vector<uint8_t> square_;
    uint8_t size_;
};

uint8_t DiagonalLatinSquare::Empty = std::numeric_limits<uint8_t>::max();

std::ostream& operator<<(std::ostream& os, const DiagonalLatinSquare& square) {
    for (uint8_t i = 0; i < square.Size(); ++i) {
        for (uint8_t j = 0; j < square.Size(); ++j) {
            os << (square(i, j) == DiagonalLatinSquare::Empty ? '.' : static_cast<char>('0' + square(i, j))) << " ";
        }
        os << std::endl;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<DiagonalLatinSquare>& squares) {
    for (size_t i = 0; i < squares.size(); ++i) {
        os << "Square " << i << ":" << std::endl;
        os << squares[i];
    }
    return os;
}

class DiagonalLatinSquareGenerationState: public DiagonalLatinSquare {
public:

    explicit DiagonalLatinSquareGenerationState(uint8_t size) : DiagonalLatinSquare(size) {
        std::vector<uint8_t> values;
        for (uint8_t i = 0; i < size; ++i) {
            values.emplace_back(i);
        }
        available_values_ = std::vector<std::vector<std::vector<uint8_t>>>(size, std::vector<std::vector<uint8_t>>(size, values));
    }

    void SetValue(uint8_t i, uint8_t j, uint8_t value) {
        if ((*this)(i, j) != Empty && (*this)(i, j) != value) {
            std::cout << "MyState " << (int)i << " " << (int)j << " " << (int)value << std::endl;
            std::cout << (*this) << std::endl;
            Application::fatalError("Trying to set different value on a non-empty position");
        }

        (*this)(i, j) = value;
        available_values_[i][j] = {value};

        for (uint8_t k = 0; k < Size(); ++k) {
            if (k != j && (*this)(i, k) == Empty) {
                RemoveFromAvailableValues(i, k, value);
                if (available_values_[i][k].empty()) {
                    promising_ &= false; return;
                }
                UpdatePossibilities(i, k);
            }
            if (k != i && (*this)(k, j) == Empty) {
                RemoveFromAvailableValues(k, j, value);
                if (available_values_[k][j].empty()) {
                    promising_ &= false; return;
                }
                UpdatePossibilities(k, j);
            }
            if (i == j && i != k && (*this)(k, k) == Empty) {
                RemoveFromAvailableValues(k, k, value);
                if (available_values_[k][k].empty()) {
                    promising_ &= false; return;
                }
                UpdatePossibilities(k, k);
            }
            if (i == Size() - j - 1 && i != k && (*this)(k, Size() - k - 1) == Empty) {
                RemoveFromAvailableValues(k, Size() - k - 1, value);
                if (available_values_[k][Size() - k - 1].empty()) {
                    promising_ &= false; return;
                }
                UpdatePossibilities(k, Size() - k - 1);
            }
        }
    }

    [[nodiscard]] bool SetAllOnes() {
        while (promising_ && !positions_with_one_possibility_.empty()) {
            auto position = std::move(positions_with_one_possibility_.back());
            positions_with_one_possibility_.pop_back();
            SetValue(position.first, position.second, available_values_[position.first][position.second].front());
        }
        return promising_;
    }

    [[nodiscard]] bool IsPromising() const {
        return promising_;
    }

    [[nodiscard]] const std::vector<uint8_t>& GetPossibleValues(uint8_t i, uint8_t j) const {
        return available_values_[i][j];
    }

    [[nodiscard]] std::vector<std::pair<uint8_t, uint8_t>> GetPositionsWithLeastPossibilities() {
        std::vector<std::pair<uint8_t, uint8_t>> positions;
        uint8_t least_possibilities = std::numeric_limits<uint8_t>::max();
        for (uint8_t i = 0; i < Size(); ++i) {
            for (uint8_t j = 0; j < Size(); ++j) {
                if ((*this)(i, j) != Empty)
                    continue;
                if (available_values_[i][j].size() < least_possibilities) {
                    positions = {{i, j}};
                    least_possibilities = available_values_[i][j].size();
                } else if (available_values_[i][j].size() == least_possibilities) {
                    positions.emplace_back(i, j);
                }
            }
        }
        return positions;
    }

    [[nodiscard]] std::vector<std::pair<uint8_t, uint8_t>> GetPositionsOnDiagonalsWithLeastPossibilities() {
        std::vector<std::pair<uint8_t, uint8_t>> positions;
        uint8_t least_possibilities = std::numeric_limits<uint8_t>::max();
        for (uint8_t i = 0; i < Size(); ++i) {
            if ((*this)(i, i) == Empty) {
                if (available_values_[i][i].size() < least_possibilities) {
                    positions = {{i, i}};
                    least_possibilities = available_values_[i][i].size();
                } else if (available_values_[i][i].size() == least_possibilities) {
                    positions.emplace_back(i, i);
                }
            }
            if ((*this)(i, Size() - i - 1) == Empty) {
                if (available_values_[i][Size() - i - 1].size() < least_possibilities) {
                    positions = {{i, Size() - i - 1}};
                    least_possibilities = available_values_[i][Size() - i - 1].size();
                } else if (available_values_[i][Size() - i - 1].size() == least_possibilities) {
                    positions.emplace_back(i, Size() - i - 1);
                }
            }
        }
        return positions;
    }

    void PrintMatrix() const {
        for (uint8_t i = 0; i < Size(); ++i) {
            for (uint8_t j = 0; j < Size(); ++j) {
                std::cout << static_cast<int>(available_values_[i][j].size()) << " ";
            }
            std::cout << std::endl;
        }
    }

    [[nodiscard]] const DiagonalLatinSquare* GetDiagonalLatinSquare() const {
        return static_cast<const DiagonalLatinSquare*>(this);
    }

    [[nodiscard]] DiagonalLatinSquare* GetDiagonalLatinSquare() {
        return static_cast<DiagonalLatinSquare*>(this);
    }

private:
    inline void RemoveFromAvailableValues(uint8_t i, uint8_t j, uint8_t value) {
        const auto& it = std::find(available_values_[i][j].begin(), available_values_[i][j].end(), value);
        if (it != available_values_[i][j].end()) {
            available_values_[i][j].erase(it);
        }
    }

    inline void RegisterOnePossibility(uint8_t i, uint8_t j) {
        auto it = std::find(positions_with_one_possibility_.begin(), positions_with_one_possibility_.end(), std::make_pair(i, j));
        if (it == positions_with_one_possibility_.end()) {
            positions_with_one_possibility_.emplace_back(i, j);
        }
    }

    inline void UpdatePossibilities(uint8_t i, uint8_t j) {
        if (available_values_[i][j].size() == 1) {
            RegisterOnePossibility(i, j);
        }
    }

private:
    std::vector<std::vector<std::vector<uint8_t> > > available_values_;
    std::vector<std::pair<uint8_t, uint8_t> > positions_with_one_possibility_;
    bool promising_ = true;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "[";
    for (auto it = v.begin(); it != v.end(); ++it) {
        os << (int)*it << " ";
    }
    os << "]";
    return os;
}

class DiagonalLatinSquareGenerator {
public:
    explicit DiagonalLatinSquareGenerator(size_t squareSize)
            : squareSize_(squareSize) {
        if ((squareSize_ > 1 && squareSize_ < 4) || squareSize_ > 15) {
            Application::fatalError("Invalid square size");
        }
    }

    [[nodiscard]] std::vector<DiagonalLatinSquare> GenerateSquares() const {
        std::vector<DiagonalLatinSquare> squares;

        std::vector<uint8_t> firstRow;
        for (uint8_t i = 0; i < squareSize_; ++i) {
            firstRow.emplace_back(i);
        }

        auto start = std::chrono::high_resolution_clock::now();
        do {
            std::vector<uint8_t> firstColumn;
            for (uint8_t i = 0; i < squareSize_; ++i)
                if (firstRow[0] != i)
                    firstColumn.emplace_back(i);
            do {
                DiagonalLatinSquareGenerationState state(squareSize_);
                for (uint8_t i = 0; i < squareSize_; ++i) {
                    state.SetValue(0, i, firstRow[i]);
                }
                bool is_valid = true;
                for (uint8_t i = 1; i < squareSize_; ++i) {
                    if (state(i, 0) != DiagonalLatinSquare::Empty && state(i, 0) != firstColumn[i - 1]) {
                        is_valid = false;
                        break;
                    }
                    state.SetValue(i, 0, firstColumn[i - 1]);
                    if (!state.SetAllOnes()) {
                        is_valid = false;
                        break;
                    }
                }
                if (!is_valid) {
                    continue;
                }
                for (auto &square: GenerateAllDiagonals(state)) {
                    if (square.SetAllOnes()) {
                        for (auto &squareState: GenerateAllPossibleSquares(std::move(square))) {
                            squares.emplace_back(std::move(squareState));
                        }
                    }
                }
                if (squares.size() > 10000000) {
                    return squares;
                }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            std::cout << static_cast<double>(squares.size()) / static_cast<double>(duration) * 1e9 << " DLS/s\tTotal Squares: " << squares.size() << std::endl;
            } while (std::next_permutation(firstColumn.begin(), firstColumn.end()));
        } while (std::next_permutation(firstRow.begin(), firstRow.end()));
        return squares;
    }

    [[nodiscard]] static std::vector<DiagonalLatinSquareGenerationState> GenerateAllDiagonals(DiagonalLatinSquareGenerationState state) {
        std::vector<DiagonalLatinSquareGenerationState> results;
        std::vector<DiagonalLatinSquareGenerationState> states;
        states.emplace_back(std::move(state));
        while (!states.empty()) {
            DiagonalLatinSquareGenerationState newState(std::move(states.back()));
            states.pop_back();
            if (!newState.IsPromising() || !newState.SetAllOnes()) {
                continue;
            }
            auto positionsOnDiagonalWithLeastPossibilities = newState.GetPositionsOnDiagonalsWithLeastPossibilities();
            if (positionsOnDiagonalWithLeastPossibilities.empty()) {
                results.emplace_back(std::move(newState));
                continue;
            }
            auto position = positionsOnDiagonalWithLeastPossibilities.front();
            for (const auto& value: newState.GetPossibleValues(position.first, position.second)) {
                DiagonalLatinSquareGenerationState currentState(newState);
                currentState.SetValue(position.first, position.second, value);
                if (currentState.IsPromising() && currentState.SetAllOnes()) {
                    states.emplace_back(std::move(currentState));
                }
            }
        }
        return results;
    }

    [[nodiscard]] static std::vector<DiagonalLatinSquareGenerationState> GenerateAllPossibleSquares(DiagonalLatinSquareGenerationState square) {
        std::vector<DiagonalLatinSquareGenerationState> results;
        std::vector<DiagonalLatinSquareGenerationState> states;
        states.emplace_back(std::move(square));
        while (!states.empty()) {
            DiagonalLatinSquareGenerationState newState(std::move(states.back()));
            states.pop_back();
            if (!newState.IsPromising() || !newState.SetAllOnes()) {
                continue;
            }
            auto positionsWithLeastPossibilities = newState.GetPositionsWithLeastPossibilities();
            if (positionsWithLeastPossibilities.empty()) {
                results.emplace_back(std::move(newState));
                continue;
            }
            auto position = positionsWithLeastPossibilities.front();
            for (const auto& value: newState.GetPossibleValues(position.first, position.second)) {
                DiagonalLatinSquareGenerationState currentState(newState);
                currentState.SetValue(position.first, position.second, value);
                if (currentState.IsPromising() && currentState.SetAllOnes()) {
                    states.emplace_back(std::move(currentState));
                }
            }
        }
        return results;
    }

private:
    uint8_t squareSize_;
};

template<typename T>
std::pair<T, T> makeOrderedPair(T a, T b) {
    if (a < b) {
        std::swap(a, b);
    }
    return std::make_pair(std::move(a), std::move(b));
}

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

class DLS_LoopFinder {
public:

    static std::pair<uint32_t, uint32_t> CountLoops(const DiagonalLatinSquare& square) {
        uint32_t fullCycles = 0;
        uint32_t partialCycles = 0;
        std::unordered_set<std::pair<uint8_t, uint8_t>, pair_hash> visited;
        for (uint8_t row = 0; row < square.Size(); ++row) {
            for (uint8_t i = 0; i < square.Size(); ++i) {
                for (uint8_t j = 0; j < square.Size(); ++j) {
                    if (i == j) continue;
                    if (visited.find(makeOrderedPair(row * square.Size() + i, row * square.Size() + j)) != visited.end()) {
                        continue;
                    }
                    uint8_t cycleLength = 1;
                    FindLoop(square, {row, i}, {row, j}, {row, i}, visited, cycleLength);
                    cycleLength == 2 * square.Size() ? ++fullCycles : ++partialCycles;
                }
            }
        }
        return {fullCycles, partialCycles};
    }

private:
    static void FindLoop(const DiagonalLatinSquare& square,
                         const std::pair<uint8_t, uint8_t>& previousPosition,
                         const std::pair<uint8_t, uint8_t>& currentPosition,
                         const std::pair<uint8_t, uint8_t>& start,
                         std::unordered_set<std::pair<uint8_t, uint8_t>, pair_hash>& visited,
                         uint8_t& length) {
        if (currentPosition == start) {
            return;
        }
        visited.insert(makeOrderedPair(currentPosition.first * square.Size() + currentPosition.second,
                                       previousPosition.first * square.Size() + previousPosition.second));
        auto value = square(previousPosition.first, previousPosition.second);
        auto position = length & 1 ? FindInColumn(square, currentPosition.second, value) : FindInRow(square, currentPosition.first, value);
        if (visited.find(makeOrderedPair(position.first * square.Size() + position.second,
                                         currentPosition.first * square.Size() + currentPosition.second)) != visited.end()) {
            return;
        }
        if (position.first == square.Size()) {
            Application::fatalError("Could not find position");
        }
        FindLoop(square, currentPosition, position, start, visited, ++length);
    }

    static std::pair<uint8_t, uint8_t> FindInRow(const DiagonalLatinSquare& square, uint8_t row, uint8_t value) {
        for (uint8_t i = 0; i < square.Size(); ++i) {
            if (square(row, i) == value) {
                return {row, i};
            }
        }
        return {square.Size(), square.Size()};
    }

    static std::pair<uint8_t, uint8_t> FindInColumn(const DiagonalLatinSquare& square, uint8_t column, uint8_t value) {
        for (uint8_t i = 0; i < square.Size(); ++i) {
            if (square(i, column) == value) {
                return {i, column};
            }
        }
        return {square.Size(), square.Size()};
    }

};

class DLS_SubrectangleFinder {
public:
    static uint64_t CountSubLatin(const DiagonalLatinSquare& square) {
        uint64_t count = 0;
        auto allPossibilities = GenerateAllSubarrays(square.Size());
        for (const auto& rowPossibility : allPossibilities) {
            for (const auto& columnPossibility: allPossibilities) {
                count += CheckSubrectangle(square, rowPossibility, columnPossibility);
            }
        }
        return count;
    }

    static bool CheckSubrectangle(const DiagonalLatinSquare& square,
                                  const std::vector<uint8_t>& rowPossibility,
                                  const std::vector<uint8_t>& columnPossibility) {

        auto size = std::max(rowPossibility.size(), columnPossibility.size());
        std::vector<uint8_t> allValues;
        // Check each row contains unique elements
        for (auto& i : rowPossibility) {
            std::vector<uint8_t> row;
            for (auto &j: columnPossibility) {
                if (std::find(row.begin(), row.end(), square(i, j)) != row.end()) {
                    return false;
                }
                if (std::find(allValues.begin(), allValues.end(), square(i, j)) == allValues.end()) {
                    allValues.push_back(square(i, j));
                    if (allValues.size() > size)
                        return false;
                }
                row.emplace_back(square(i, j));
            }
        }

        // Check each column contains unique elements
        for (auto& j : columnPossibility) {
            std::vector<uint8_t> column;
            for (auto& i : rowPossibility) {
                if (std::find(column.begin(), column.end(), square(i, j)) != column.end()) {
                    return false;
                }
                column.emplace_back(square(i, j));
            }
        }

        return true;
    }

    static std::vector<std::vector<uint8_t>> GenerateAllSubarrays(uint8_t n, bool generateEmpty = false) {
        std::vector<std::vector<uint8_t>> result;
        std::queue<std::pair<std::vector<uint8_t>, uint8_t>> queue;
        queue.push({{}, 0});
        while (!queue.empty()) {
            auto current = std::move(queue.front());
            queue.pop();
            if (current.second == n) {
                if (generateEmpty ||!current.first.empty()) {
                    result.push_back(current.first);
                }
            } else {
                queue.push({current.first, current.second + 1});
                current.first.emplace_back(current.second++);
                queue.emplace(std::move(current));
            }
        }
        return result;
    }
};

class SequenceGenerator {
public:
    template<typename T>
    static std::vector<T> GenerateSequence(T n, uint64_t k) {
        // Fast O(n) solution
        uint64_t factorial = 1;
        std::vector<T> nums;
        for (T i = 0; i < n; ++i) {
            nums.emplace_back(i);
            factorial *= (i + 1);
        }
        std::vector<T> result;
        for (T i = 0; i < n - 1; ++i) {
            factorial /= (n - i);
            T current = k / factorial;
            k = k % factorial;
            result.emplace_back(std::move(nums[current]));
            nums.erase(nums.begin() + current);
        }
        result.emplace_back(std::move(nums[0]));
        return result;
    }

    template<typename T>
    static std::vector<T> GenerateSequence(std::vector<std::vector<T>> possibleValues, uint64_t k) {
        // Slow but more flexible solution
        std::queue<std::vector<T>> queue;
        queue.push({});
        while (!queue.empty()) {
            auto current = std::move(queue.front());
            queue.pop();
            if (current.size() == possibleValues.size()) {
                if (k == 0) {
                    return current;
                }
                --k;
            } else {
                for (auto& value : possibleValues[current.size()]) {
                    if (std::find(current.begin(), current.end(), value) == current.end()) {
                        auto newCurrent = current;
                        newCurrent.emplace_back(value);
                        queue.push(std::move(newCurrent));
                    }
                }
            }
        }
        throw -1;
        return {};
    }

    template<typename T>
    static size_t CountSequences(std::vector<std::vector<T>> possibleValues) {
        // Slow but more flexible solution
        size_t count = 0;
        std::queue<std::vector<T>> queue;
        queue.push({});
        while (!queue.empty()) {
            auto current = std::move(queue.front());
            queue.pop();
            if (current.size() == possibleValues.size()) {
                ++count;
            } else {
                for (auto& value : possibleValues[current.size()]) {
                    if (std::find(current.begin(), current.end(), value) == current.end()) {
                        auto newCurrent = current;
                        newCurrent.emplace_back(value);
                        queue.push(std::move(newCurrent));
                    }
                }
            }
        }
        return count;
    }
};
