/* =============================================================================

 Copyright (C) 2009-2021 Valerii Sukhorukov. All Rights Reserved.

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

================================================================================
*/
/**
* \file in_loop_action.h
* Eables run-time updates of control loop parameters.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_IN_LOOP_ACTION_H
#define UTILS_IN_LOOP_ACTION_H

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <type_traits>
#include <utility>

#include "msgr.h"

/// Library outer namespace.
namespace utils {

/**
 * \class InLoopAction
 * \brief Implements updating of loop control parameters in iterative process.
 * \tparam Q Type of the parameter variable.
 */
template <typename Q>
class InLoopAction {

    template <typename S, typename Enabler=void>
    class Field {};

    template <typename S>
    class Field<S*, typename std::enable_if_t<std::is_pointer_v<S*>>> {

    public:

        const char* name;

        Field(const std::pair<S*, const char*>& f)
            : name {f.second}
            , val {f.first}
        {}

        auto& operator()() noexcept { return val; }

        void initialize( std::ofstream& ofs ) const
        {
            if (val)
                ofs << (val ? *val : 0)
                    << "\t\t\t\t\t# " << name
                    << std::endl;
        }

        void read_new( std::ifstream& fin )
        {
            if (!val) return;

            std::string line;
            std::getline(fin, line);
            std::stringstream ss {line};
            ss >> newval;
        }

        auto adjust( const std::string& pName )
        {
            if (!val) return false;

            if (newval != *val) {    // take the new field over
                std::cout << "\nApplying new value " << name << ": "<< newval
                          << " to the "+pName
                          << std::endl
                          << "old one is: " << *val
                          << std::endl;

                *val = newval;
                return true;
            }
            return false;
        }

        void print()
        {
            std::cout << "\t" << name << " "
                      << (val ? *val : 0)
                      << std::endl;
        }

    private:

        S* val {};
        S newval {};
    };

    template <typename S>
    class Field<S, typename std::enable_if_t<!std::is_pointer_v<S>>> {

    public:

        const S defval;
        const char* name;

        Field( S defval,
               const char* name )
            : defval {defval}
            , name {name}
        {}

        auto& operator()() noexcept { return val; }

        void initialize( std::ofstream& ofs ) const
        {
            ofs << val
                << std::string("\t\t\t\t\t# ")+name
                << std::endl;
        }

        void read_new( std::ifstream& fin )
        {
            std::string line;
            std::getline(fin, line);
            std::stringstream ss {line};
            ss >> newval;
        }

        auto adjust( const std::string& pName )
        {
            if (newval != val) {    // take the new field over
                std::cout << "\nApplying new " << name << " "<< newval
                          << " to the "+pName
                          << std::endl
                          << "old one is: " << val << std::endl;

                val = newval;
                return true;
            }
            return false;
        }

        void print()
        {
            std::cout << "\t" << name << " "
                      << val
                      << std::endl;
        }

    private:

        S val {defval};
        S newval {};
    };

public:

    Field<Q> start;
    Field<Q> end;
    Field<Q> delta;
    Field<Q> next {huge<Q>, "next"};

    Field<bool*> flag0;
    Field<bool*> flag1;
    Field<bool*> flag2;

    explicit InLoopAction ( Msgr& msgr,
                const std::filesystem::path& path,
                const std::string& s,
                const std::string& parName,
                const Q& val,        // const ref
                const Q startDefault,
                const Q endDefault,
                const Q deltaDefault,
                const std::pair<bool*, const char*>& flag0={},
                const std::pair<bool*, const char*>& flag1={},
                const std::pair<bool*, const char*>& flag2={}
                );

    void check_file();
    void update_next() noexcept;

private:

    const std::filesystem::path path;
    const std::filesystem::path file;

    const std::string parName;

    volatile const Q& val;

    Msgr& msgr;

    void initialize_checkFile() const;
};

template <typename Q>
InLoopAction<Q>::
InLoopAction ( Msgr& msgr,
                const std::filesystem::path& path,
                const std::string& s,
                const std::string& parName,
                const Q& val,        // const ref
                const Q startDefault,
                const Q endDefault,
                const Q deltaDefault,
                const std::pair<bool*, const char*>& flag0,
                const std::pair<bool*, const char*>& flag1,
                const std::pair<bool*, const char*>& flag2
                 )
    : start {startDefault, "start"}
    , end {endDefault, "end"}
    , delta {deltaDefault, "delta"}
    , flag0 {flag0}
    , flag1 {flag1}
    , flag2 {flag2}
    , path {path}
    , file {path / s}
    , parName {parName}
    , val {val}
    , msgr {msgr}
{
    if (!std::filesystem::exists(path))
        std::filesystem::create_directory(path);
    initialize_checkFile();
    update_next();
}

template <typename Q>
void InLoopAction<Q>::
initialize_checkFile() const
{
    std::ofstream ofs {file};
    if (ofs.fail())
        msgr.exit("Unable to create check file ", file);

    start.initialize(ofs);
    end.initialize(ofs);
    delta.initialize(ofs);
    flag0.initialize(ofs);
    flag1.initialize(ofs);
    flag2.initialize(ofs);
}

template <typename Q>
void InLoopAction<Q>::
check_file()
{
    std::ifstream fin {file};
    if (fin.fail())
        msgr.exit("Unable to open file for reading: ", file);

    start.read_new(fin);
    end.read_new(fin);
    delta.read_new(fin);
    flag0.read_new(fin);
    flag1.read_new(fin);
    flag2.read_new(fin);

    fin.close();

    const auto a = start.adjust(parName);
    const auto b = end.adjust(parName);
    const auto c = delta.adjust(parName);
    flag0.adjust(parName);
    flag1.adjust(parName);
    flag2.adjust(parName);

    if (a || b || c) update_next();        // flags are adjusted inside and does not affect _next
}

template <typename Q>
void InLoopAction<Q>::
update_next() noexcept
{
    std::cout << "In "+parName << std::endl;
    start.print();
    end.print();
    delta.print();

    next() = (start() >= val)
           ? (start() <= end())
             ? start()
             : huge<Q>
          : (val + delta() > end())
             ? huge<Q>
             : val + delta();
    std::cout << "Next event is at: " << next() << std::endl
               << std::endl;
}

}  // namespace utils

#endif  // UTILS_IN_LOOP_ACTION_H
