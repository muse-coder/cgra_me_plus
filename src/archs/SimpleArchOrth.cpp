#include <CGRA/Exception.h>
#include <CGRA/user-inc/UserArchs.h>
#include <CGRA/user-inc/UserModules.h>
std::unique_ptr<CGRA> UserArchs::createSimpleArchOrth(const std::map<std::string, std::string> & args)
{
    int cols;
    int rows;
    int homo;
    try
    {
        cols = std::stoi(args.at("cols"));
        rows = std::stoi(args.at("rows"));
        homo = std::stoi(args.at("homogeneous_fu"));
    }
    catch(const std::out_of_range & e)
    {
        throw cgrame_error("C++ Architecture Argument Error");
    }

    std::unique_ptr<CGRA> result(new CGRA());

    // create top IOs
    for (int i = 0; i < cols; i++)
    {
        result->addPort("ext_io_top_" + std::to_string(i), PORT_BIDIR);
        result->addSubModule(new IO("io_top_" + std::to_string(i), 32));
        result->addConnection(
                "io_top_" + std::to_string(i) + ".bidir",
                "this.ext_io_top_" + std::to_string(i)
        ); // to top-level external pin
    }

    // create right IOs
    for (int i = 0; i < rows; i++)
    {
        result->addPort("ext_io_right_" + std::to_string(i), PORT_BIDIR);
        result->addSubModule(new IO("io_right_" + std::to_string(i), 32));
        result->addConnection(
                "io_right_" + std::to_string(i) + ".bidir",
                "this.ext_io_right_" + std::to_string(i)
        ); // to top-level external pin
    }

    // create bottom IOs
    for (int i = 0; i < cols; i++)
    {
        result->addPort("ext_io_bottom_" + std::to_string(i), PORT_BIDIR);
        result->addSubModule(new IO("io_bottom_" + std::to_string(i), 32));
        result->addConnection(
                "io_bottom_" + std::to_string(i) + ".bidir",
                "this.ext_io_bottom_" + std::to_string(i)
        ); // to top-level external pin
    }

    // create left IOs
    for (int i = 0; i < rows; i++)
    {
        result->addPort("ext_io_left_" + std::to_string(i), PORT_BIDIR);
        result->addSubModule(new IO("io_left_" + std::to_string(i), 32));
        result->addConnection(
                "io_left_" + std::to_string(i) + ".bidir",
                "this.ext_io_left_" + std::to_string(i)
        ); // to top-level external pin
    }

    // create mem column
    for (int i = 0; i < rows; i++)
    {
        result->addSubModule(new MemPort("mem_" + std::to_string(i), cols)); // create enough input ports for the outputs of each FU in a row
    }


    // create cols and rows of FUs
    if(homo)
    {
        for (int c = 0; c < cols; c++)
        {
            for (int r = 0; r < rows; r++)
            {
                result->addSubModule(new SimpleArchFU("b_c" + std::to_string(c) + "_r" + std::to_string(r), 5));
            }
        }
    }
    else // hetero checker pattern
    {
        for (int c = 0; c < cols; c++)
        {
            for (int r = 0; r < rows; r++)
            {
                if((c + r) % 2) // create
                {
                    result->addSubModule(new SimpleArchFU("b_c" + std::to_string(c) + "_r" + std::to_string(r), 5, false));
                }
                else
                {
                    result->addSubModule(new SimpleArchFU("b_c" + std::to_string(c) + "_r" + std::to_string(r), 5, true));
                }
            }
        }
    }

    // create all connections

    // top row I/O connections
    for(int c = 0; c < cols; c++)
    {
        std::string io_name = "io_top_" + std::to_string(c);
        std::string block_name = "b_c" + std::to_string(c) + "_r" + std::to_string(0);
        result->addConnection(io_name + ".out", block_name + ".in0"); // to block
        result->addConnection(block_name + ".out", io_name + ".in"); // to io
    }

    // right col I/O connections
    for(int r = 0; r < rows; r++)
    {
        std::string io_name = "io_right_" + std::to_string(r);
        std::string block_name = "b_c" + std::to_string(cols-1) + "_r" + std::to_string(r);
        result->addConnection(io_name + ".out", block_name + ".in1"); // to block
        result->addConnection(block_name + ".out", io_name + ".in"); // to io
    }

    // bottom row I/O connections
    for(int c = 0; c < cols; c++)
    {
        std::string io_name = "io_bottom_" + std::to_string(c);
        std::string block_name = "b_c" + std::to_string(c) + "_r" + std::to_string(rows-1);
        result->addConnection(io_name + ".out", block_name + ".in2"); // to block
        result->addConnection(block_name + ".out", io_name + ".in"); // to io
    }

    // left col I/O connections
    for(int r = 0; r < rows; r++)
    {
        std::string io_name = "io_left_" + std::to_string(r);
        std::string block_name = "b_c" + std::to_string(0) + "_r" + std::to_string(r);
        result->addConnection(io_name + ".out", block_name + ".in3"); // to block
        result->addConnection(block_name + ".out", io_name + ".in"); // to io
    }

    // MEM connections
    // full crossbar for each row
    for(int r = 0; r < rows; r++)
    {
        std::string mem_name = "mem_" + std::to_string(r);

        for(int c = 0; c < cols; c++)
        {
            std::string block_name = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
            result->addConnection(block_name + ".out", mem_name + ".in" + std::to_string(c)); // to io
            result->addConnection(mem_name + ".out", block_name + ".in4"); // to block
        }
    }

    // internal north/south connections
    for (int c = 0; c < cols; c++)
    {
        for (int r = 0; r < rows - 1; r++)
        {
            std::string block_name_c_r  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
            std::string block_name_c_r1 = "b_c" + std::to_string(c) + "_r" + std::to_string(r + 1);

            // north / south connections
            result->addConnection(block_name_c_r + ".out", block_name_c_r1 + ".in0");
            result->addConnection(block_name_c_r1 + ".out", block_name_c_r + ".in2");
        }
    }

    // internal west/east connections
    for (int c = 0; c < cols - 1; c++)
    {
        for (int r = 0; r < rows; r++)
        {
            std::string block_name_c_r  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
            std::string block_name_c1_r = "b_c" + std::to_string(c + 1) + "_r" + std::to_string(r);

            // east / west connections
            result->addConnection(block_name_c_r + ".out", block_name_c1_r + ".in3");
            result->addConnection(block_name_c1_r + ".out", block_name_c_r + ".in1");
        }
    }

    return result;
}
