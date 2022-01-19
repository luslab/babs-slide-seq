/*
 * nourdine.bah@crick.ac.uk
 */

#include <string>
#include <vector>
#include "types.hpp"

std::vector<std::string> get_sequences(char* path);
Sequence seq_to_num(std::string sequence);
std::string num_to_seq(Sequence number, int length);
char distance(Sequence number1, Sequence number2, char length);
int get_length(std::vector<std::string> reads, std::vector<std::string> barcodes, char* length);
