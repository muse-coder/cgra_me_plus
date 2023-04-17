/*******************************************************************************
 * CGRA-ME Software End-User License Agreement
 *
 * The software programs comprising "CGRA-ME" and the documentation provided
 * with them are copyright by its authors S. Chin, K. Niu, N. Sakamoto, J. Zhao,
 * A. Rui, S. Yin, A. Mertens, J. Anderson, and the University of Toronto. Users
 * agree to not redistribute the software, in source or binary form, to other
 * persons or other institutions. Users may modify or use the source code for
 * other non-commercial, not-for-profit research endeavours, provided that all
 * copyright attribution on the source code is retained, and the original or
 * modified source code is not redistributed, in whole or in part, or included
 * in or with any commercial product, except by written agreement with the
 * authors, and full and complete attribution for use of the code is given in
 * any resulting publications.
 *
 * Only non-commercial, not-for-profit use of this software is permitted. No
 * part of this software may be incorporated into a commercial product without
 * the written consent of the authors. The software may not be used for the
 * design of a commercial electronic product without the written consent of the
 * authors. The use of this software to assist in the development of new
 * commercial CGRA architectures or commercial soft processor architectures is
 * also prohibited without the written consent of the authors.
 *
 * This software is provided "as is" with no warranties or guarantees of
 * support.
 *
 * This Agreement shall be governed by the laws of Province of Ontario, Canada.
 *
 * Please contact Prof. Anderson if you are interested in commercial use of the
 * CGRA-ME framework.
 ******************************************************************************/

#include <string>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <memory>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <chrono>

#include <linux/limits.h>
#include <unistd.h>
#include <libgen.h>

#include <CGRA/user-inc/UserArchs.h>

#include <CGRA/Exception.h>
#include <CGRA/Mapper.h>
#include <CGRA/Mapping.h>

#include <CGRA/dotparse.h>
#include <CGRA/adlparse.h>

#include <CGRA/dotparse.h>

#include <CGRA/Visual.h>

#include <cxxopts.hpp>

#include <mini.hpp>

#include <CGRA/lisa/gnn.h>
#include <CGRA/lisa/LISAController.h>
#include <CGRA/lisa/LISAMapper.h>

#include <sys/types.h>
#include <sys/stat.h>

//FIXME: actually, we do not know which edge is backedge in CGRA_ME, though we can find the cycle.
    // But cgra-me seems order the node that somehow follow the data dependency. Let us try to use this now.


void find_backedge(std::vector<std::pair<int, int>> lisa_edges, std::set<std::pair<int, int>> lisa_backedges, std::map<int, OpGraphOp*> & id_to_node)
{
    std::map<int, std::set<int>> node_transitive_children; //
    std::map<int, std::set<int>> node_transitive_parents;  //
    std::set<int> global_visited_node;

    for (auto info : id_to_node)
    {
        node_transitive_children.emplace(info.first, std::set<int>());
        node_transitive_parents.emplace(info.first, std::set<int>());
    }


    //FIXME: this cannot detect all the backedges. Please compare with the one in LISAController:init labels.
    for (auto e : lisa_edges)
    {
        if (node_transitive_children[e.second].find(e.first) != node_transitive_children[e.second].end())
        {
            //backedge
            lisa_backedges.insert(std::make_pair(e.first, e.second));
        }
        else
        {
            // std::cout<<"visit"<<e.src<<","<<e.des<<std::endl;
            for (auto parent : node_transitive_parents[e.first])
            {
                node_transitive_parents[e.second].insert(parent);
            }
            node_transitive_children[e.first].insert(e.second);
            node_transitive_parents[e.second].insert(e.first);
            for (auto child : node_transitive_children[e.second])
            {
                node_transitive_children[e.first].insert(child);
            }
        }
    }
    std::vector<std::pair<int, int>> temp_edges;
    for(auto e: lisa_edges){
        if(lisa_backedges.find(e) == lisa_backedges.end()){
            temp_edges.push_back(e);
        }
    }
    assert(temp_edges.size() + lisa_backedges.size() == lisa_edges.size());
    lisa_edges = temp_edges;
}

void set_edges(std::shared_ptr<OpGraph> opgraph_, std::vector<std::pair<int, int>> &lisa_edges, std::set<std::pair<int, int>> & lisa_backedges,  std::map<OpGraphOp*, int>  node_to_id_,  std::map<int, OpGraphOp*> & id_to_node)
{

    LOG(LISAGNN) << "dump edges in original operator graph";

    // this is from cgra_me
    // std::map<OpGraphOp *, int> opnode_map;
    // int counter = 0;
    // for (auto it = opgraph_->op_nodes.begin(); it != opgraph_->op_nodes.end(); it++)
    // {
    //     assert(*it);
    //     opnode_map[(*it)] = counter++;
    //     // std::cout << opnode_map[(*it)] << "[opcode=" << (*it)->opcode << "];\n";
    // }

    for (auto it = opgraph_->val_nodes.begin(); it != opgraph_->val_nodes.end(); it++)
    {
        int inputnode = node_to_id_[(*it)->input];

        assert((*it)->output.size() == (*it)->output_operand.size());
        for (unsigned int o = 0; o < (*it)->output.size(); o++)
        {
            OpGraphOp *op = (*it)->output[o];
            unsigned int operand = (*it)->output_operand[o];
            if(inputnode != node_to_id_[op] )lisa_edges.push_back(std::make_pair(inputnode, node_to_id_[op]));
            // std::cout << inputnode << "->" << opnode_map[op] << "[operand=" << operand << "]; ";
            // std::cout << "//" << (*it)->input->name << "->" << op->name << "\n";
        }
    }

    LOG(LISAGNN) << " edges";
    for (auto it : lisa_edges)
    {
        LOG(LISAGNN) << "  (" << it.first << "," << it.second << ")";
    }

    ;
    find_backedge(lisa_edges, lisa_backedges, id_to_node);
    LOG(LISAGNN) << "back edges";
    for (auto it : lisa_backedges)
    {
        LOG(LISAGNN) << "  (" << it.first << "," << it.second << ")";
    }
}

int main(int argc, char* argv[])
{
    std::cout << std::endl;
    std::cout << "CGRA - Modelling and Exploration Version 1.0 (http://cgra-me.ece.utoronto.ca/)" << std::endl;
    std::cout << "Copyright (c) 2015-2018 University of Toronto. All Rights Reserved." << std::endl;
    std::cout << "For research and academic purposes only. Commercial use is prohibited." << std::endl;
    std::cout << "Please email questions to: Xander Chin(xan@ece.utoronto.ca)" << std::endl;
    std::cout << "Compiled: " << __DATE__ << " " << __TIME__ << std::endl;
    std::cout << std::endl;

     for(int i=0; i < argc; i++){
        std::cout << "Argument " << i << ": " << argv[i] << std::endl;
    }

    UserArchs userarchs;

    std::string arch_filename, dfg_filename, hdl_dirpath;
    int arch_id;
    std::string arch_opts;
    MapperType mapper_type;
    std::string mapper_opts;
    int max_II;
    int start_II;
    double timelimit;
    bool printarch;
    bool printop;
    bool genverilog;
    bool make_testbench;
    int adl;
    int cgra_x;
    int cgra_y;
    bool lisa_training;
    bool lisa_inef;
    bool lisa_eval_routing_priority;
    bool lisa_eval_unmapped;
    std::string training_output_file;
    std::string lisa_dfg_id = "none"; // this is for LISA and dump final labels
    std::string arch_model_name = "cgra_me"; //used as folder name of dumped label
    bool parallel_rouing_and_place = false;
    bool SA_more_running_time = false;

    try
    {
        cxxopts::Options options("CGRA-ME", "CGRA - Modelling and Exploration");

        options.add_options()
            ("h,help", "Print Help")
            ("x,xml", "Use Architecture Description Language File", cxxopts::value<std::string>(), "<Filepath>")
            ("p,parser", "Which ADL parser to use (0 = XML specs v0, 1 = XML specs v1)", cxxopts::value<int>()->default_value("1"), "<#>")
            ("c,cpp", "Use C++ Architecture, ID # Generated from --arch-list", cxxopts::value<int>(), "<#>")
            ("arch-list", "Show the List of Avaliable C++ Architectures with IDs", cxxopts::value<bool>())
            ("arch-opts", "C++ Architecture Options that Overwrites the Default Ones (<Key>=<Value> Pairs, Separate by Space, and Close by Quotation Marks)", cxxopts::value<std::string>(), "<\"opts\">")
            ("arch-opts-list", "Show the List of Avaliable Options for a C++ Architecture, ID # Generated from --arch-list", cxxopts::value<int>(), "<#>")
            ("i,II", "Architecture Contexts, Zhaoying: this should be max II value", cxxopts::value<int>()->default_value("32"), "<#>")
            ("sII", "start II", cxxopts::value<int>()->default_value("1"), "<#>")
            ("g,dfg", "The DFG file to map in dot format", cxxopts::value<std::string>())
            ("dfg_id", "The is for LISA training, and dump ", cxxopts::value<std::string>()->default_value("none"))
            ("m,mapper", "Which Mapper to Use (0 = ILP, 1 = Simulated Annealing, 2 = LISA mapper)", cxxopts::value<int>()->default_value("0"), "<#>")
            ("mapper-opts", "Mapper Options that Overwrites the Default Ones (<Key>=<Value> Pairs, Separate by Space, and Close by Quotation Marks)", cxxopts::value<std::string>(), "<\"opts\">")
            ("v,visual", "Output visualization directory after mapping", cxxopts::value<bool>())
            ("t,timelimit", "Mapper Timelimit", cxxopts::value<double>()->default_value("7200.0"), "<#>")
            ("training", "training for LISA", cxxopts::value<bool>())
            ("inef", "inference with LISA", cxxopts::value<bool>())
            ("eval_routing_priority", "evaluate routing priority", cxxopts::value<bool>())
            ("eval_unmapped", "evaluate unmapped", cxxopts::value<bool>())
            ("ppr", "parallel routing and placement", cxxopts::value<bool>())
            ("arch_name", "arch name",  cxxopts::value<std::string>()->default_value("cgra_me"), "<#>")
            ("training_output", "training output log for LISA",  cxxopts::value<std::string>()->default_value("lisa_training"), "<#>")
            ("a,print-arch", "Print Architecture to stdout", cxxopts::value<bool>())
            ("o,print-op", "Print Operation Graph to stdout", cxxopts::value<bool>())
            ("sa_more_running_time", "add more SA running time", cxxopts::value<bool>())
            ("cgra_x", "cgra_x", cxxopts::value<int>()->default_value("12"), "<#>")
            ("cgra_y", "cgra_y", cxxopts::value<int>()->default_value("4"), "<#>")
            ("gen-verilog", "Generate Verilog Implementation of Architecture and Dump to Specified Directory", cxxopts::value<std::string>(), "<Directorypath>")
            ("gen-testbench", "Generate testbench for use in a simulation of the DFG with configuration bitstream", cxxopts::value<bool>())
            ;

        options.parse(argc, argv);

        genverilog = options.count("gen-verilog");

        if(options.count("help")) // Print help
        {
            std::cout << options.help({""}) << std::endl;
            return 0;
        }

        if(options.count("arch-list")) // Print the list of avaliable CGRA architecture
        {
            std::cout << "List of C++ CGRA Architecture with ID: " << std::endl;
            int count = 0;
            for(const auto & arch : userarchs.archs)
                std::cout << count++ << ". " << arch.second.first << std::endl;
            return 0;
        }

        if(options.count("arch-opts-list")) // Print the list of avaliable options for a CGRA architecture
        {
            int arch_index(options["arch-opts-list"].as<int>());
            if(!(arch_index >= 0 && arch_index < userarchs.archs.size()))
            {
                std::cout << "[ERROR] The C++ Architecture ID Given Is out of Bound, Run with --arch-list to See the List of Archs with IDs" << std::endl;
                return 1;
            }
            const auto & arch_entry = *std::next(userarchs.archs.begin(), arch_index);
            const auto & arch_default_args = arch_entry.second.second;
            std::cout << "List of Options for: " << arch_entry.second.first << std::endl;
            std::cout << "Option:" << std::setw(80) << "Default Value:" << std::endl;
            for(const auto & arch_arg : arch_default_args)
                std::cout << std::left << std::setw(40) << arch_arg.first << std::right << std::setw(40) << arch_arg.second << std::endl;
            return 1;
        }

        if(!genverilog && !options.count("dfg"))
        {
            //std::cout << "[ERROR] DFG File Path is Missing, Exiting..." << std::endl;
            std::cout << options.help({""}) << std::endl;
            return 1;
        }

        if(options.count("xml") && options.count("cpp"))
        {
            std::cout << "[ERROR] You Can Only Choose Either ADL or C++ but Not Both" << std::endl;
            return 1;
        }
        else if(!options.count("xml") && !options.count("cpp"))
        {
            std::cout << "[ERROR] No Architecture Specified, You Can Use Either ADL or C++" << std::endl;
            return 1;
        }

        if(options.count("cpp"))
        {
            int cpp_index(options["cpp"].as<int>());
            if(!(cpp_index >= 0 && cpp_index < userarchs.archs.size()))
            {
                std::cout << "[ERROR] The C++ Architecture ID Given Is out of Bound, Run with --arch-list to See the List of Archs with IDs" << std::endl;
                return 1;
            }
        }

        dfg_filename = options["dfg"].as<std::string>();
        arch_filename = options["xml"].as<std::string>();
        hdl_dirpath = options["gen-verilog"].as<std::string>();
        arch_id = options["cpp"].as<int>();
        arch_opts = options["arch-opts"].as<std::string>();
        mapper_type = static_cast<MapperType>(options["mapper"].as<int>());
        mapper_opts = options["mapper-opts"].as<std::string>();
        max_II = options["II"].as<int>();
        start_II = options["sII"].as<int>();
        timelimit = options["timelimit"].as<double>();
        printarch = options["print-arch"].as<bool>();
        printop = options["print-op"].as<bool>();
        adl = options["parser"].as<int>();
        make_testbench = options["gen-testbench"].as<bool>();
        cgra_x = options["cgra_x"].as<int>();
        cgra_y = options["cgra_y"].as<int>();
        lisa_training = options["training"].as<bool>();
        lisa_inef = options["inef"].as<bool>();
        lisa_eval_unmapped = options["eval_unmapped"].as<bool>();
        lisa_eval_routing_priority = options["eval_routing_priority"].as<bool>();
        training_output_file = options["training_output"].as<std::string>();
        arch_model_name = options["arch_name"].as<std::string>();
        lisa_dfg_id = options["dfg_id"].as<std::string>();
        parallel_rouing_and_place = options["ppr"].as<bool>();
        SA_more_running_time = options["sa_more_running_time"].as<bool>();  
    }
    catch(const cxxopts::OptionException & e)
    {
        std::cout << "[ERROR] Error Parsing Options: " << e.what() << std::endl;
        return 1;
    }

    try
    {
        char full_path[PATH_MAX] = {0}; // initialize to zero, so we don't have to set the null char ourselves
        ssize_t count = readlink("/proc/self/exe", full_path, PATH_MAX-1);
        std::string exe_path;
        if(count != -1)
        {
            exe_path = std::string(dirname(full_path));
        }
        else
        {
            std::cout << "[ERROR] Readlink is not able to get the executable path" << std::endl;
            return 1;
        }

        std::shared_ptr<CGRA> arch = NULL;

        if(arch_filename.empty()) // Use C++
        {
            std::cout << "[INFO] Creating Architecture #" << arch_id << " from C++..." << std::endl;
            auto & arch_entry = *std::next(userarchs.archs.begin(), arch_id);
            std::cout << "[INFO] Architecture Name: " << arch_entry.second.first << std::endl;
            auto & arch_default_args = arch_entry.second.second;
            std::stringstream ss_arch_opts(arch_opts);
            std::string opt_pair;
            while(ss_arch_opts >> opt_pair)
            {
                auto n = opt_pair.find('=');
                if(n == std::string::npos)
                {
                    std::cout <<"[WARNING] Ill Formated C++ Architrecture Option, Skipping: " << opt_pair << std::endl;
                    continue;
                }
                auto key = opt_pair.substr(0, n);
                auto value = opt_pair.substr(n + 1, std::string::npos);
                auto it = arch_default_args.find(key);
                if(it != arch_default_args.end())
                {
                    std::cout << "[INFO] Overwritting C++ Architecture Parameter: " << key << " to " << value << " (Default: " << arch_default_args.at(key) << ")" << std::endl;
                    arch_default_args[key] = value; // Overwrite default option
                }
                else
                    std::cout << "[WARNING] C++ Architecture Parameter: " << key << " doesn't exist, Skipping: " << key << " = " << value << std::endl;
            }
            std::cout << "[INFO] Creating \"" << arch_entry.second.first << "\" Architecture from C++..." << std::endl;
            arch = arch_entry.first(arch_default_args); // Create C++ Arch
        }
        else // Use ADL
        {
            std::cout << "[INFO] Creating Architecture from XML..." << std::endl;
            if(adl)
            {
                std::string adl_template_filename = exe_path + "/module_templates.xml";
                std::ifstream adl_template_file(adl_template_filename);
                if (adl_template_file.good())
                    arch = adl1::parseADL(adl_template_filename, arch_filename);
                else
                {
                    std::cout << "[ERROR] Missing module_templates.xml, Cannot Read Default Module Templates" << std::endl;
                    return 1;
                }
            }
            else
                arch = adl0::parseADL(arch_filename);
        }

        // Print Architecture if requested
        if(printarch)
        {
            std::cout << std::endl << "[ARCHGRAPH]" << std::endl;
            arch->print();
            arch->print_dot();
        }

        // Generate Verilog implementation and exit
        if(genverilog)
        {
            // Verify if directory exists
            struct stat fs_info;
            auto res = stat(hdl_dirpath.c_str(), &fs_info);
            if (res != 0)
            {
                std::cout << "[ERROR] Cannot open \"" + hdl_dirpath + "\"";
                std::cout << std::endl;
                return 1;
            }
            bool is_dir = fs_info.st_mode & S_IFDIR;
            if (!is_dir) // make sure it is directory
            {
                std::cout << "[ERROR] \"" + hdl_dirpath + "\" is not a directory";
                std::cout << std::endl;
                return 1;
            }

            std::cout << "[INFO] Generating Verilog Implementation of Specified Architecture" << std::endl;
            if (hdl_dirpath.back() != '/')
                hdl_dirpath = hdl_dirpath + "/";
            arch->genVerilog(hdl_dirpath);
            return 0;
        }

        // Creating OpGraph
        std::cout << "[INFO] Parsing DFG..." << std::endl;
        // Need to do std::move for GCC 6.2 (possibly others)
        std::shared_ptr<OpGraph> opgraph = parseOpGraph(dfg_filename.c_str());


        // Print OpGraph if requested
        if(printop)
        {
            std::cout << std::endl << "[OPGRAPH]" << std::endl;
            opgraph->print_dot();
        }

        std::ifstream ini_file(exe_path + "/mapper_config.ini");
        if(!ini_file)
        {
            std::cout << "[ERROR] Missing mapper_config.ini, Cannot Read Default Parameters for Mappers" << std::endl;
            return 1;
        }
        std::string ini_str((std::istreambuf_iterator<char>(ini_file)), std::istreambuf_iterator<char>());

        mINI config_ini;
        std::map<std::string, std::string> mapper_args;
        if(config_ini.parse(ini_str)) // INI Parse OK
        {
            mapper_args = config_ini;
        }
        else // INI Parse ERROR
        {
            std::cout <<"[ERROR] Error in INI Parser, Check You mapper_config.ini" << std::endl;
            return 1;
        }

        std::stringstream ss_mapper_opts(mapper_opts);
        std::string opt_pair;
        while(ss_mapper_opts >> opt_pair)
        {
            auto n = opt_pair.find('=');
            if(n == std::string::npos)
            {
                std::cout <<"[WARNING] Ill Formated Mapper Option, Skipping: " << opt_pair << std::endl;
                continue;
            }
            auto key = opt_pair.substr(0, n);
            auto value = opt_pair.substr(n + 1, std::string::npos);
            auto it = mapper_args.find(key);
            if(it != mapper_args.end())
            {
                std::cout << "[INFO] Overwritting Mapper Parameter: " << key << " to " << value << " (Default: " << mapper_args.at(key) << ")" << std::endl;
                mapper_args[key] = value; // Overwrite default option
            }
            else
                std::cout << "[WARNING] Mapper Parameter: " << key << " doesn't exist, Skipping: " << key << " = " << value << std::endl;
        }
        if(SA_more_running_time){
            // assert(mapper_args.find("AnnealMapper.swap_factor")!=mapper_args.end());
            // std::cout << "[INFO] allow SA to search more steps...  old:"<<mapper_args["AnnealMapper.swap_factor"];
            // mapper_args["AnnealMapper.swap_factor"] = std::to_string(std::stoi( mapper_args["AnnealMapper.swap_factor"]) * 10 );
            // std::cout<<" new:"<<mapper_args["AnnealMapper.swap_factor"]<<std::endl;


            // assert(mapper_args.find("AnnealMapper.cold_accept_rate")!=mapper_args.end());
            // std::cout << "[INFO] allow explore all temperature...  old:"<<mapper_args["AnnealMapper.cold_accept_rate"];
            // mapper_args["AnnealMapper.cold_accept_rate"] = std::to_string(0);
            // std::cout<<" new:"<<mapper_args["AnnealMapper.cold_accept_rate"]<<std::endl;
            timelimit = 3600 * 4;

        }
        std::cout << "[INFO] Creating Mapper..." << std::endl;
        arch->ROWS = cgra_x;
        arch->COLS = cgra_y;
        
        

        

        
        // simply calculated by DFG/fu node
        MRRG* mrrg = arch->getMRRG(1).get();
        for (auto &fu : mrrg->function_nodes)
        {
            std::cout<<fu->getFullName()<<"\n";
        }
        // assert(false);
        int fu_num = (cgra_x -2)* (cgra_y -2);// cgra-me did not provide an elegant way to get the fu/pe number;
        
        std::set<OpGraphOp*> unique_set;
        for(auto node:opgraph->op_nodes ){
            if(node->opcode !=OpGraphOpCode::OPGRAPH_OP_INPUT && node->opcode !=OpGraphOpCode::OPGRAPH_OP_OUTPUT){
                unique_set.insert(node);
            }
        }
        int dfg_num  = unique_set.size();
        int min_II = ceil( (float)dfg_num / fu_num);
        std::cout<<"fu num:"<<fu_num<< " dfg num:"<<dfg_num<< " minII:"<<min_II<<"\n";
        std::string result_filename = "result/";
        if(lisa_eval_routing_priority){
            result_filename = result_filename + "eval_routing_priority_"+ arch_model_name+ "_result.txt";
        }else if(lisa_eval_unmapped){
            result_filename = result_filename + "eval_unmapped_"+ arch_model_name+ "_result.txt";
        }else if(SA_more_running_time){
            result_filename = result_filename + "sa_more_running_time_"+ arch_model_name+ "_result.txt";
        }else{
            result_filename = result_filename + arch_model_name+ "_result.txt";
        }
        
        std::ofstream result_file;
        

        std::string method;
        if(mapper_type == MapperType::ILPMapper){
            method = "ILP";
        }else if(mapper_type ==  MapperType::AnnealMapper){
             method = "SA";
        }else if(mapper_type == MapperType::LISAMapper){
            method = "i-LISA";
            if(parallel_rouing_and_place){
                method = "para_i-LISA"; // no need to use this one
            }
        }


            

         

        int running_time = 1;
        if(mapper_type == MapperType::AnnealMapper || lisa_eval_routing_priority ){
            running_time = 3;
        }
        std::vector<perf_metric> perf_hist;
        for(int i = 0; i < running_time; i++ ){

            int currII = start_II;
            currII = currII < min_II? min_II: currII;
            
            auto start = std::chrono::steady_clock::now();

    
            while(currII <= max_II){
                std::cout<<"*************** try II:"<<currII<<" max II:"<<max_II<<std::endl;
                auto mapper = Mapper::createMapper(mapper_type, arch, timelimit, mapper_args);

    
                Mapping mapping_result = mapper->mapOpGraph(opgraph, currII, arch_model_name);

                if(mapping_result.isMapped())
                {
                    std::cout << "*************** mapped successful on II "<<currII<<" *********************"<<std::endl;
                    genMappingVisual(exe_path, mapping_result);
                    mapping_result.outputMapping();

                    if (make_testbench)
                    {
                        std::ofstream tb_file("testbench.v");
                        arch->genBitStream(mapping_result).print_testbench(tb_file);
                    }

                    break;
                }
                std::cout << "*************** mapped failed on II "<<currII<<std::endl;
                currII ++;
            }
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;
            int used_second =   elapsed_seconds.count() ;
            std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
            perf_metric this_iter_perf { currII, 0 , used_second };
            perf_hist.push_back(this_iter_perf);
        }
        std::sort(perf_hist.begin(), perf_hist.end());

        int median_index = 0;
        for(int j = 0; j < running_time; j++){
            std::cout<<"perf"<<j <<" "<<perf_hist[j].ii<<""<<perf_hist[j].running_time<<"\n";
        }
        if(running_time == 3) {
            median_index = 1;
            assert(perf_hist.size() == 3);
        }
        if(SA_more_running_time){
            assert(mapper_type ==  MapperType::AnnealMapper );
            method = "SA_T";
        }
        perf_metric final_perf = perf_hist[median_index];
        result_file.open (result_filename, std::ios_base::app); 
        result_file << arch_filename<<" "<<dfg_filename<<" "<<method<<" "<<min_II<<" "<<final_perf.ii<<" "<<final_perf.running_time<<"\n";
        result_file.close();
        
        return 0; // cannot map
    }
    catch(const cgrame_error & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 1;
};

