#pragma once
#include "../support/types.hpp"

namespace nustar::star2711b::solo::detail {

class ParameterBase {
public:
    std::string nameString; // string that identifies parameter
    int inputLevel;         // where the parameter was defined
    int inputLevelAllowed;  // at which input level parameter definition is allowed
    virtual void inputValues(std::istringstream& streamIn) = 0;
    friend std::ostream& operator<<(std::ostream& o, ParameterBase const& b);
    virtual ~ParameterBase() {};

protected:
    virtual void printValues(std::ostream& o) const = 0;
};

inline std::ostream& operator<<(std::ostream& o, ParameterBase const& b) {
    b.printValues(o);
    return o;
};

template <class parameterType> inline parameterType inputOneValue(std::istringstream& streamIn) {
    parameterType oneV;
    streamIn >> oneV;
    return oneV;
};
template <> inline std::string inputOneValue<std::string>(std::istringstream& streamIn) {
    std::string oneV = "";
    streamIn >> std::ws;          // skip whitespace
    if (streamIn.peek() != '"') { // simple parameter with no spaces or "
        streamIn >> oneV;
    } else {
        streamIn.get(); // skip "
        std::getline(streamIn, oneV, '"');
    };
    return oneV;
};

template <class parameterType>
inline void printOneValue(parameterType* value, std::ostream& outStr) {
    outStr << *value;
};
template <> inline void printOneValue<std::string>(std::string* value, std::ostream& outStr) {
    if ((*value).find_first_of(" \t") !=
        std::string::npos) { // there is white space in the argument, put "" around
        outStr << '\"' << *value << '\"';
    } else {
        outStr << *value;
    };
};

template <class parameterType> class ScalarParameter : public ParameterBase {
public:
    parameterType* value;
    std::vector<parameterType> allowedValues;

    ScalarParameter(int inputLevelIn, int inputLevelAllowedIn, std::string nameStringIn,
                    parameterType* valueIn) {
        nameString = nameStringIn;
        inputLevel = inputLevelIn;
        inputLevelAllowed = inputLevelAllowedIn;
        value = valueIn;
    };

    void inputValues(std::istringstream& streamIn) {
        *value = inputOneValue<parameterType>(streamIn);
    };

    ~ScalarParameter() {};

protected:
    virtual void printValues(std::ostream& outStr) const {
        printOneValue(value, outStr);
    };
};

template <class parameterType> class VectorParameter : public ParameterBase {
public:
    std::vector<parameterType>* value;
    std::vector<parameterType> allowedValues;

    VectorParameter(int inputLevelIn, int inputLevelAllowedIn, std::string nameStringIn,
                    std::vector<parameterType>* valueIn) {
        nameString = nameStringIn;
        inputLevel = inputLevelIn;
        inputLevelAllowed = inputLevelAllowedIn;
        value = valueIn;
    };

    void inputValues(std::istringstream& streamIn) {
        (*value).clear();
        while (streamIn.good()) {
            (*value).push_back(inputOneValue<parameterType>(streamIn));
            streamIn >> std::ws; // remove white space, may arrive at the end of line
        };
    };

    ~VectorParameter() {};

protected:
    virtual void printValues(std::ostream& outStr) const {
        for (int ii = 0; ii < (int)(*value).size(); ii++) {
            printOneValue(&(*value).at(ii), outStr);
            outStr << "   ";
        };
    };
};

} // namespace nustar::star2711b::solo::detail
