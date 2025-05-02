#include <iostream>

#include "argparse.h"
#include "library.h"
#include "abm.h"


int main(int argc, char* argv[]) {
    argparse::ArgumentParser main_program("abm");

    main_program.add_description("Agent Based Modelling");

    main_program.add_argument("--edgelist")
        .required()
        .help("Input edgelist (source, target)");
    main_program.add_argument("--nodelist")
        .required()
        .help("Input nodelist (node, year)");
    main_program.add_argument("--out-degree-bag")
        .required()
        .help("Input out-degree bag (, out-degree)");
    main_program.add_argument("--recency-probabilities")
        .required()
        .help("Input recency bag (year, probability)");
    main_program.add_argument("--alpha")
        .default_value(double(-1))
        .help("Neighborhood alpha")
        .scan<'g', double>();
    main_program.add_argument("--growth-rate")
        .required()
        .help("Growth rate")
        .scan<'g', double>();
    main_program.add_argument("--num-cycles")
        .required()
        .help("Number of years")
        .scan<'d', int>();
    main_program.add_argument("--output-file")
        .required()
        .help("Output clustering file");
    main_program.add_argument("--log-file")
        .required()
        .help("Output log file");
    main_program.add_argument("--num-processors")
        .default_value(int(1))
        .help("Number of processors")
        .scan<'d', int>();
    main_program.add_argument("--log-level")
        .default_value(int(1))
        .help("Log level where 0 = silent, 1 = info, 2 = verbose")
        .scan<'d', int>();
    try {
        main_program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << main_program;
        std::exit(1);
    }

    std::string edgelist = main_program.get<std::string>("--edgelist");
    std::string nodelist = main_program.get<std::string>("--nodelist");
    std::string out_degree_bag = main_program.get<std::string>("--out-degree-bag");
    std::string recency_probabilities = main_program.get<std::string>("--recency-probabilities");
    double alpha = main_program.get<double>("--alpha");
    double growth_rate = main_program.get<double>("--growth-rate");
    int num_cycles = main_program.get<int>("--num-cycles");
    std::string output_file = main_program.get<std::string>("--output-file");
    std::string log_file = main_program.get<std::string>("--log-file");
    int num_processors = main_program.get<int>("--num-processors");
    int log_level = main_program.get<int>("--log-level") - 1; // so that enum is cleaner
    ABM* abm = new ABM(edgelist, nodelist, out_degree_bag, recency_probabilities, alpha, growth_rate, num_cycles, output_file, log_file, num_processors, log_level);
    abm->main();
    delete abm;
}
