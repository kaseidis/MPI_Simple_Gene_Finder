/**
 * @file InputParser.h
 * @author iain@stackoverflow
 * @brief Argument parser.
 *        Reference: https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
 *
 * @copyright iain@stackoverflow
 */
#pragma once
#ifndef _INPUT_PARSER_H
#define _INPUT_PARSER_H
#include <string>
#include <vector>

class InputParser
{
public:
    InputParser(int &argc, char **argv);
    /// @author iain
    const std::string &getCmdOption(const std::string &option) const;
    /// @author iain
    bool cmdOptionExists(const std::string &option) const;

private:
    std::vector<std::string> tokens;
};
#endif