//#include <CGRA/Exception.h>
//#include <CGRA/user-inc/UserArchs.h>
//#include <CGRA/user-inc/UserModules.h>
//
//std::unique_ptr<CGRA> UserArchs::createSimpleArchDiag(const std::map<std::string, std::string> & args)
//{
//    int cols;
//    int rows;
//    int homo;
//    try
//    {
//        cols = std::stoi(args.at("cols"));
//        rows = std::stoi(args.at("rows"));
//        homo = std::stoi(args.at("homogeneous_fu"));
//    }
//    catch(const std::out_of_range & e)
//    {
//        throw cgrame_error("C++ Architecture Argument Error");
//    }
//
//    std::unique_ptr<CGRA> result(new CGRA());
//
//    // create corner IOs
//    for (int i = 0; i < 4; i++)
//    {
//        result->addPort("ext_io_corner_" + std::to_string(i), PORT_BIDIR);
//        result->addSubModule(new IOPort("io_corner_" + std::to_string(i), 1));
//        result->addConnection(
//                "io_corner_" + std::to_string(i) + ".bidir",
//                "this.ext_io_corner_" + std::to_string(i)
//        ); // to top-level external pin
//    }
//
//    // create top IOs
//    for (int i = 0; i < cols; i++)
//    {
//        result->addPort("ext_io_top_" + std::to_string(i), PORT_BIDIR);
//        if(i == 0 || i == cols-1)
//        {
//            result->addSubModule(new IOPort("io_top_" + std::to_string(i), 2));
//        }
//        else
//        {
//            result->addSubModule(new IOPort("io_top_" + std::to_string(i), 3));
//        }
//        result->addConnection(
//                "io_top_" + std::to_string(i) + ".bidir",
//                "this.ext_io_top_" + std::to_string(i)
//        ); // to top-level external pin
//    }
//
//    // create right IOs
//    for (int i = 0; i < rows; i++)
//    {
//        result->addPort("ext_io_right_" + std::to_string(i), PORT_BIDIR);
//        if(i == 0 || i == rows-1)
//        {
//            result->addSubModule(new IOPort("io_right_" + std::to_string(i), 2));
//        }
//        else
//        {
//            result->addSubModule(new IOPort("io_right_" + std::to_string(i), 3));
//        }
//        result->addConnection(
//                "io_right_" + std::to_string(i) + ".bidir",
//                "this.ext_io_right_" + std::to_string(i)
//        ); // to top-level external pin
//    }
//
//    // create bottom IOs
//    for (int i = 0; i < cols; i++)
//    {
//        result->addPort("ext_io_bottom_" + std::to_string(i), PORT_BIDIR);
//        if(i == 0 || i == cols-1)
//        {
//            result->addSubModule(new IOPort("io_bottom_" + std::to_string(i), 2));
//        }
//        else
//        {
//            result->addSubModule(new IOPort("io_bottom_" + std::to_string(i), 3));
//        }
//        result->addConnection(
//                "io_bottom_" + std::to_string(i) + ".bidir",
//                "this.ext_io_bottom_" + std::to_string(i)
//        ); // to top-level external pin
//    }
//
//    // create left IOs
//    for (int i = 0; i < rows; i++)
//    {
//        result->addPort("ext_io_left_" + std::to_string(i), PORT_BIDIR);
//        if(i == 0 || i == rows-1)
//        {
//            result->addSubModule(new IOPort("io_left_" + std::to_string(i), 2));
//        }
//        else
//        {
//            result->addSubModule(new IOPort("io_left_" + std::to_string(i), 3));
//        }
//        result->addConnection(
//                "io_left_" + std::to_string(i) + ".bidir",
//                "this.ext_io_left_" + std::to_string(i)
//        ); // to top-level external pin
//    }
//
//    // create mem ports
//    for (int i = 0; i < rows; i++)
//    {
//        result->addSubModule(new MemPort("mem_" + std::to_string(i), cols)); // create enough input ports for the outputs of each FU in a row
//    }
//
//    // create cols and rows of FUs
//    if(homo)
//    {
//        for (int c = 0; c < cols; c++)
//        {
//            for (int r = 0; r < rows; r++)
//            {
//                result->addSubModule(new SimpleArchFU("b_c" + std::to_string(c) + "_r" + std::to_string(r), 9));
//            }
//        }
//    }
//    else // hetero checker pattern
//    {
//        for (int c = 0; c < cols; c++)
//        {
//            for (int r = 0; r < rows; r++)
//            {
//                if((c + r) % 2) // create
//                {
//                    result->addSubModule(new SimpleArchFU("b_c" + std::to_string(c) + "_r" + std::to_string(r), 9, false));
//                }
//                else
//                {
//                    result->addSubModule(new SimpleArchFU("b_c" + std::to_string(c) + "_r" + std::to_string(r), 9, true));
//                }
//            }
//        }
//
//    }
//
//    // create all connections
//
//    // top row I/O connections for each I/0
//    for(int c = 0; c < cols; c++)
//    {
//        std::string io_name = "io_top_" + std::to_string(c);
//        std::string block_name = "b_c" + std::to_string(c) + "_r" + std::to_string(0);
//        std::string block_name_p1 = "b_c" + std::to_string(c+1) + "_r" + std::to_string(0);
//        std::string block_name_m1 = "b_c" + std::to_string(c-1) + "_r" + std::to_string(0);
//
//        if(c == 0)
//        {
//            result->addConnection(io_name + ".out", block_name + ".in2"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name_p1 + ".in1"); // to block
//            result->addConnection(block_name_p1 + ".out", io_name + ".in1"); // to io
//        }
//        else if(c == cols-1)
//        {
//            result->addConnection(io_name + ".out", block_name_m1 + ".in3"); // to block
//            result->addConnection(block_name_m1 + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name + ".in2"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in1"); // to io
//        }
//        else
//        {
//            result->addConnection(io_name + ".out", block_name_m1 + ".in3"); // to block
//            result->addConnection(block_name_m1 + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name + ".in2"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in1"); // to io
//
//            result->addConnection(io_name + ".out", block_name_p1 + ".in1"); // to block
//            result->addConnection(block_name_p1 + ".out", io_name + ".in2"); // to io
//        }
//    }
//
//    // right col I/O connections
//    for(int r = 0; r < rows; r++)
//    {
//        std::string io_name = "io_right_" + std::to_string(r);
//        std::string block_name = "b_c" + std::to_string(cols-1) + "_r" + std::to_string(r);
//        std::string block_name_m1 = "b_c" + std::to_string(cols-1) + "_r" + std::to_string(r-1);
//        std::string block_name_p1 = "b_c" + std::to_string(cols-1) + "_r" + std::to_string(r+1);
//
//        if(r == 0)
//        {
//            result->addConnection(io_name + ".out", block_name + ".in4"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name_p1 + ".in3"); // to block
//            result->addConnection(block_name_p1 + ".out", io_name + ".in1"); // to io
//        }
//        else if(r == rows-1)
//        {
//            result->addConnection(io_name + ".out", block_name_m1 + ".in5"); // to block
//            result->addConnection(block_name_m1 + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name + ".in4"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in1"); // to io
//        }
//        else
//        {
//            result->addConnection(io_name + ".out", block_name_m1 + ".in5"); // to block
//            result->addConnection(block_name_m1 + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name + ".in4"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in1"); // to io
//
//            result->addConnection(io_name + ".out", block_name_p1 + ".in3"); // to block
//            result->addConnection(block_name_p1 + ".out", io_name + ".in2"); // to io
//        }
//    }
//
//    // bottom row I/O connections
//    for(int c = 0; c < cols; c++)
//    {
//        std::string io_name = "io_bottom_" + std::to_string(c);
//        std::string block_name = "b_c" + std::to_string(c) + "_r" + std::to_string(rows-1);
//        std::string block_name_p1 = "b_c" + std::to_string(c+1) + "_r" + std::to_string(rows-1);
//        std::string block_name_m1 = "b_c" + std::to_string(c-1) + "_r" + std::to_string(rows-1);
//
//        if(c == 0)
//        {
//            result->addConnection(io_name + ".out", block_name + ".in6"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name_p1 + ".in7"); // to block
//            result->addConnection(block_name_p1 + ".out", io_name + ".in1"); // to io
//        }
//        else if(c == cols-1)
//        {
//            result->addConnection(io_name + ".out", block_name_m1 + ".in5"); // to block
//            result->addConnection(block_name_m1 + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name + ".in6"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in1"); // to io
//        }
//        else
//        {
//            result->addConnection(io_name + ".out", block_name_m1 + ".in5"); // to block
//            result->addConnection(block_name_m1 + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name + ".in6"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in1"); // to io
//
//            result->addConnection(io_name + ".out", block_name_p1 + ".in7"); // to block
//            result->addConnection(block_name_p1 + ".out", io_name + ".in2"); // to io
//        }
//    }
//
//    // left col I/O connections
//    for(int r = 0; r < rows; r++)
//    {
//        std::string io_name = "io_left_" + std::to_string(r);
//        std::string block_name = "b_c" + std::to_string(0) + "_r" + std::to_string(r);
//        std::string block_name_m1 = "b_c" + std::to_string(0) + "_r" + std::to_string(r-1);
//        std::string block_name_p1 = "b_c" + std::to_string(0) + "_r" + std::to_string(r+1);
//
//        if(r == 0)
//        {
//            result->addConnection(io_name + ".out", block_name + ".in0"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name_p1 + ".in1"); // to block
//            result->addConnection(block_name_p1 + ".out", io_name + ".in1"); // to io
//        }
//        else if(r == rows-1)
//        {
//            result->addConnection(io_name + ".out", block_name_m1 + ".in7"); // to block
//            result->addConnection(block_name_m1 + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name + ".in0"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in1"); // to io
//        }
//        else
//        {
//            result->addConnection(io_name + ".out", block_name_m1 + ".in7"); // to block
//            result->addConnection(block_name_m1 + ".out", io_name + ".in0"); // to io
//
//            result->addConnection(io_name + ".out", block_name + ".in0"); // to block
//            result->addConnection(block_name + ".out", io_name + ".in1"); // to io
//
//            result->addConnection(io_name + ".out", block_name_p1 + ".in1"); // to block
//            result->addConnection(block_name_p1 + ".out", io_name + ".in2"); // to io
//        }
//    }
//
//    // corner I/O connections
//    {
//        std::string io_name;
//        std::string block_name;
//        // NW
//        io_name = "io_corner_" + std::to_string(0);
//        block_name = "b_c" + std::to_string(0) + "_r" + std::to_string(0);
//        result->addConnection(io_name + ".out", block_name + ".in1"); // to block
//        result->addConnection(block_name + ".out", io_name + ".in0"); // to io
//        // NE
//        io_name = "io_corner_" + std::to_string(1);
//        block_name = "b_c" + std::to_string(cols-1) + "_r" + std::to_string(0);
//        result->addConnection(io_name + ".out", block_name + ".in3"); // to block
//        result->addConnection(block_name + ".out", io_name + ".in0"); // to io
//        // SW
//        io_name = "io_corner_" + std::to_string(2);
//        block_name = "b_c" + std::to_string(cols-1) + "_r" + std::to_string(rows-1);
//        result->addConnection(io_name + ".out", block_name + ".in5"); // to block
//        result->addConnection(block_name + ".out", io_name + ".in0"); // to io
//        // SE
//        io_name = "io_corner_" + std::to_string(3);
//        block_name = "b_c" + std::to_string(0) + "_r" + std::to_string(rows-1);
//        result->addConnection(io_name + ".out", block_name + ".in7"); // to block
//        result->addConnection(block_name + ".out", io_name + ".in0"); // to io
//    }
//
//    // MEM connections
//    // full crossbar for each row
//    for(int r = 0; r < rows; r++)
//    {
//        std::string mem_name = "mem_" + std::to_string(r);
//
//        for(int c = 0; c < cols; c++)
//        {
//            std::string block_name = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
//            result->addConnection(block_name + ".out", mem_name + ".in" + std::to_string(c)); // to io
//            result->addConnection(mem_name + ".out", block_name + ".in8"); // to block
//        }
//    }
//
//    // internal perimeter connections
//    // iterate over perimiter blocks and make all outbound connections to other blocks (not IOs) for each block
//
//    // top row to W, E, SE, S, SW
//    for (int c = 1; c < cols-1; c++)
//    {
//        int r = 0;
//
//        std::string block_name  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
//
//        std::string block_name_W = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r);
//        result->addConnection(block_name + ".out", block_name_W + ".in4"); // to block to the W
//
//        std::string block_name_E = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r);
//        result->addConnection(block_name + ".out", block_name_E + ".in0"); // to block to the E
//
//        std::string block_name_SE = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_SE + ".in1"); // to block to the SE
//
//        std::string block_name_S = "b_c" + std::to_string(c) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_S + ".in2"); // to block to the S
//
//        std::string block_name_SW = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_SW + ".in3"); // to block to the SW
//    }
//
//    // bottom row to W, E, NW, N, NE
//    for (int c = 1; c < cols-1; c++)
//    {
//        int r = rows-1;
//
//        std::string block_name  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
//
//        std::string block_name_W = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r);
//        result->addConnection(block_name + ".out", block_name_W + ".in4"); // to block to the W
//
//        std::string block_name_E = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r);
//        result->addConnection(block_name + ".out", block_name_E + ".in0"); // to block to the E
//
//        std::string block_name_NW = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_NW + ".in5"); // to block to the NW
//
//        std::string block_name_N = "b_c" + std::to_string(c) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_N + ".in6"); // to block to the N
//
//        std::string block_name_NE = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_NE + ".in7"); // to block to the NE
//    }
//
//    // left col to N, S, NE, E, SE
//    for (int r = 1; r < rows-1; r++)
//    {
//        int c = 0;
//
//        std::string block_name  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
//
//        std::string block_name_N = "b_c" + std::to_string(c) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_N + ".in6"); // to block to the N
//
//        std::string block_name_S = "b_c" + std::to_string(c) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_S + ".in2"); // to block to the S
//
//        std::string block_name_NE = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_NE + ".in7"); // to block to the NE
//
//        std::string block_name_E = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r);
//        result->addConnection(block_name + ".out", block_name_E + ".in0"); // to block to the E
//
//        std::string block_name_SE = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_SE + ".in1"); // to block to the SE
//    }
//
//    // right col to N, S, NW, W, SW
//    for (int r = 1; r < rows-1; r++)
//    {
//        int c = cols-1;
//
//        std::string block_name  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
//
//        std::string block_name_N = "b_c" + std::to_string(c) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_N + ".in6"); // to block to the N
//
//        std::string block_name_S = "b_c" + std::to_string(c) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_S + ".in2"); // to block to the S
//
//        std::string block_name_W = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r);
//        result->addConnection(block_name + ".out", block_name_W + ".in4"); // to block to the W
//
//        std::string block_name_NW = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_NW + ".in5"); // to block to the NW
//
//        std::string block_name_SW = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_SW + ".in3"); // to block to the SW
//    }
//
//    // corners
//    // NW to E, S, SE
//    {
//        int c = 0;
//        int r = 0;
//
//        std::string block_name  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
//
//        std::string block_name_S = "b_c" + std::to_string(c) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_S + ".in2"); // to block to the S
//
//        std::string block_name_E = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r);
//        result->addConnection(block_name + ".out", block_name_E + ".in0"); // to block to the E
//
//        std::string block_name_SE = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_SE + ".in1"); // to block to the SE
//    }
//
//    // NE to W, S, SW
//    {
//        int c = cols-1;
//        int r = 0;
//
//        std::string block_name  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
//
//        std::string block_name_S = "b_c" + std::to_string(c) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_S + ".in2"); // to block to the S
//
//        std::string block_name_W = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r);
//        result->addConnection(block_name + ".out", block_name_W + ".in4"); // to block to the W
//
//        std::string block_name_SW = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r+1);
//        result->addConnection(block_name + ".out", block_name_SW + ".in3"); // to block to the SW
//    }
//
//    // SW to N, E, NE
//    {
//        int c = 0;
//        int r = rows-1;
//
//        std::string block_name  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
//
//        std::string block_name_N = "b_c" + std::to_string(c) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_N + ".in6"); // to block to the N
//
//        std::string block_name_NE = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_NE + ".in7"); // to block to the NE
//
//        std::string block_name_E = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r);
//        result->addConnection(block_name + ".out", block_name_E + ".in0"); // to block to the E
//    }
//
//    // SE to N, W, NW
//    {
//        int c = cols-1;
//        int r = rows-1;
//
//        std::string block_name  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
//
//        std::string block_name_N = "b_c" + std::to_string(c) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_N + ".in6"); // to block to the N
//
//        std::string block_name_W = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r);
//        result->addConnection(block_name + ".out", block_name_W + ".in4"); // to block to the W
//
//        std::string block_name_NW = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r-1);
//        result->addConnection(block_name + ".out", block_name_NW + ".in5"); // to block to the NW
//    }
//
//    // iterate over each internal block making all outbound connections
//    for (int c = 1; c < cols-1; c++)
//    {
//        for (int r = 1; r < rows-1; r++)
//        {
//            std::string block_name  = "b_c" + std::to_string(c) + "_r" + std::to_string(r);
//
//            std::string block_name_W = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r);
//            result->addConnection(block_name + ".out", block_name_W + ".in4"); // to block to the W
//
//            std::string block_name_NW = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r-1);
//            result->addConnection(block_name + ".out", block_name_NW + ".in5"); // to block to the NW
//
//            std::string block_name_N = "b_c" + std::to_string(c) + "_r" + std::to_string(r-1);
//            result->addConnection(block_name + ".out", block_name_N + ".in6"); // to block to the N
//
//            std::string block_name_NE = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r-1);
//            result->addConnection(block_name + ".out", block_name_NE + ".in7"); // to block to the NE
//
//            std::string block_name_E = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r);
//            result->addConnection(block_name + ".out", block_name_E + ".in0"); // to block to the E
//
//            std::string block_name_SE = "b_c" + std::to_string(c+1) + "_r" + std::to_string(r+1);
//            result->addConnection(block_name + ".out", block_name_SE + ".in1"); // to block to the SE
//
//            std::string block_name_S = "b_c" + std::to_string(c) + "_r" + std::to_string(r+1);
//            result->addConnection(block_name + ".out", block_name_S + ".in2"); // to block to the S
//
//            std::string block_name_SW = "b_c" + std::to_string(c-1) + "_r" + std::to_string(r+1);
//            result->addConnection(block_name + ".out", block_name_SW + ".in3"); // to block to the SW
//
//        }
//    }
//
//    return result;
//}
