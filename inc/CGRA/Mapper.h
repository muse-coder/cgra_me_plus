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

#ifndef ___MAPPER_H__
#define ___MAPPER_H__

#include <map>
#include <vector>
#include <iostream>

#include <CGRA/CGRA.h>
#include <CGRA/OpGraph.h>
#include <CGRA/debug.h>
#include <CGRA/Mapping.h>


#define systolic_array_x 7
#define systolic_array_y 7
struct mpos{
    int x;
    int y;
    unsigned int t;
    bool operator==(const mpos& other) const {
    return x == other.x && y == other.y && t == other.t;
    }

     bool operator<(const mpos& other) const {
    return x < other.x || y < other.y ;
    }
};


enum class MapperType
{
    ILPMapper = 0,
    AnnealMapper = 1,
    LISAMapper = 2
};

typedef struct
{
    std::map<OpGraphNode*, std::vector<MRRGNode*>> mapping;
} OpMapping;

class Mapper
{
    public:
        virtual Mapping mapOpGraph(std::shared_ptr<OpGraph> opgraph);
        virtual Mapping mapOpGraph(std::shared_ptr<OpGraph> opgraph, int II, std::string arch_model_name) = 0;
        void genBitstream();



        virtual ~Mapper();
        static std::unique_ptr<Mapper> createMapper(MapperType mt, std::shared_ptr<CGRA> cgra, int timelimit, const std::map<std::string, std::string> & args);

    protected:
        Mapper(std::shared_ptr<CGRA> cgra, int timelimit);
        std::string arch_model_name_;

        std::map<std::pair<int, int>, int> systolic_pe_index;
        
        int sys_arr[systolic_array_x][systolic_array_y] =
        {
            {0,0,0,0,0,0,0},
            {0,1,2,1,2,2,0},
            {0,2,1,1,1,1,0},
            {0,2,1,2,1,1,0},
            {0,2,1,1,2,1,0},
            {0,2,1,2,1,1,0},
            {0,0,0,0,0,0,0}
        };

        //  int sys_arr[systolic_array_x][systolic_array_y] =
        // {
        //     {0,0,0,0,0,0,0, 0},
        //     {0,1,2,2,2,2,2, 0},
        //     {0,2,1,1,1,1,1, 0},
        //     {0,2,1,1,1,1,1, 0},
        //     {0,2,1,1,2,1,2, 0},
        //     {0,2,1,1,1,1,2, 0},
        //     {0,2,2,1,2,1,2, 0},
        //     {0,0,0,0,0,0,1, 0}
        // };

        std::shared_ptr<CGRA>   cgra;       // Architecture Object
        int     timelimit;  // The mapper timeout in seconds
};

#endif

