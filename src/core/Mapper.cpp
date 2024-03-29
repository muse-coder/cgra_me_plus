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

#include <CGRA/CGRA.h>

#include <CGRA/Mapper.h>
#include <CGRA/Mapping.h>
#include <CGRA/ILPMapper.h>
#include <CGRA/AnnealMapper.h>
#include <CGRA/lisa/LISAMapper.h>

#include <CGRA/Exception.h>

std::unique_ptr<Mapper> Mapper::createMapper(MapperType mt, std::shared_ptr<CGRA> cgra, int timelimit, const std::map<std::string, std::string> & args)
{
   
    
    switch(mt)
    {
        case MapperType::ILPMapper:
            return std::make_unique<ILPMapper>(cgra, timelimit, args);
        case MapperType::AnnealMapper:
            return std::make_unique<AnnealMapper>(cgra, timelimit, args);
        case MapperType::LISAMapper:
            return std::make_unique<LISAMapper>(cgra, timelimit, args);
    }
    throw cgrame_error("Invalid Mapper Type Specified");
}

Mapper::Mapper(std::shared_ptr<CGRA> cgra, int timelimit)
{
    this->cgra = cgra;
    this->timelimit = timelimit;
}

Mapper::~Mapper()
{
}

Mapping Mapper::mapOpGraph(std::shared_ptr<OpGraph> opgraph)
{
    assert(false && "should not call this function");
    int ii = 1;
    bool mapped = false;
    Mapping result(cgra, ii, opgraph);
    while(ii <= cgra->maxII)
    {
        result  =  mapOpGraph(opgraph, ii, "");
        if(result.isMapped())
            return result;        
    }

    return result;
}

bool operator<(const std::pair<float,MRRGNode*> & a, const std::pair<float,MRRGNode*> & b)
{
    return a.first > b.first;
}
