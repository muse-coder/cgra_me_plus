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

#include <list>
#include <algorithm>
#include <utility>
#include <queue>
#include <functional>

#include <assert.h>
#include <sys/time.h>

#include <CGRA/Exception.h>
#include <CGRA/CGRA.h>
#include <CGRA/OpGraph.h>
#include <CGRA/lisa/LISAMapper.h>

using std::cout;
using std::endl;

#define ALLOW_MULTIPLE_PLACEMENT
#define SIMPLE_ANNEAL_SCHEDULE
#define CALCULATE_INITIAL_TEMPERATURE
// Debugging output flags
//#define ANNEAL_DEBUG
//#define PLACEMENT_DEBUG
// #define DEBUG_ROUTING

LISAMapper::LISAMapper(std::shared_ptr<CGRA> cgra, int timelimit, const std::map<std::string, std::string> &args)
    : Mapper(cgra, timelimit)
{
    try
    {
        rand_seed = std::stoi(args.at("AnnealMapper.random_seed"));
        pfactor = std::stof(args.at("AnnealMapper.initial_pfactor"));
        pfactor_factor = std::stof(args.at("AnnealMapper.pfactor_factor"));
        const_temp_factor = std::stof(args.at("AnnealMapper.constant_temp_factor"));
        swap_factor = std::stoi(args.at("AnnealMapper.swap_factor"));
        cold_accept_rate = std::stof(args.at("AnnealMapper.cold_accept_rate"));
    }
    catch (const std::exception &e)
    {
        throw cgrame_error(std::string("AnnealMapper Parameter Parsing Exception Thrown by: [") + e.what() + "] at File: " + std::string(__FILE__) + " Line: " + std::to_string(__LINE__));
    }
}
LISAMapper::LISAMapper(std::shared_ptr<CGRA> cgra, int timelimit, int rand_seed, float initial_pfactor, float pfactor_factor, float const_temp_factor, int swap_factor, float cold_accept_rate)
    : Mapper(cgra, timelimit)
{
    this->rand_seed = rand_seed;
    // Initial Penalty Factor
    this->pfactor = initial_pfactor; //0.001;
    this->pfactor_factor = pfactor_factor;
    this->const_temp_factor = const_temp_factor;
    this->swap_factor = swap_factor;
    this->cold_accept_rate = cold_accept_rate;
}

LISAMapper::LISAMapper(std::shared_ptr<CGRA> cgra)
    : LISAMapper(cgra, 0, 0, 0.001, 1.05, 0.99, 100, 0.01)
{
}

float LISAMapper::getCost(MRRGNode *n)
{
    float base_cost;

    if (n->type == MRRG_NODE_FUNCTION)
        base_cost = 2.0;
    else
        base_cost = 1.0;

    return base_cost * occupancy[n] + (occupancy[n] <= n->capacity ? 0.0 : (occupancy[n] - n->capacity) * pfactor);
}

float LISAMapper::getCost(MRRG *mrrg)
{
    float total = 0.0;
    for (auto &node : mrrg->function_nodes)
    {
        total += getCost(node);
    }

    for (auto &node : mrrg->routing_nodes)
    {
        total += getCost(node);
    }

    return total;
}

float LISAMapper::getFinalCost()
{
    float total = 0.0;
    for (auto &node : mrrg_->function_nodes)
    {
        total += getCost(node);
    }

    for (auto &node : mrrg_->routing_nodes)
    {
        total += getCost(node);
    }

    return total;
}

float LISAMapper::getCost(OpGraph *opgraph)
{
    float result = 0.0;
    for (auto &op : opgraph->op_nodes)
    {
        result += getCost(op);
    }
    for (auto &val : opgraph->val_nodes)
    {
        result += getCost(val);
    }
    return result;
}

float LISAMapper::getCost(OpGraphNode *n)
{
    float result = 0.0;

    if (mapping[n].size() == 0)
    {
        cout << "there's an unroute/unmap for OpGraphNode: " << n->name << "\n";
        result = INFINITY;
    }
    else
    {
        for (auto &node : mapping[n])
        {
            result += getCost(node);
        }
    }
    return result;
}

bool LISAMapper::checkOveruse(MRRG *mrrg)
{
    bool result = true;
    overuse_fu.clear();
    for (auto &node : mrrg->function_nodes)
    {
        if (occupancy[node] > node->capacity)
        {
            LOG(OVERUSE) << *node << " fu is overused. (" << occupancy[node] << "/" << node->capacity << ")\n";
            result = false;
            overuse_fu.insert(node);
        }
    }

    for (auto &node : mrrg->routing_nodes)
    {
        if (occupancy[node] > node->capacity)
        {
            
            LOG(OVERUSE) << *node << " is overused. (" << occupancy[node] << "/" << node->capacity << ")\n";
            // for(auto mappinged_node: occupancy_detail[node]){
            //     LOG(OVERUSE) <<"use "<<mappinged_node->name;
            // }
            //  LOG(OVERUSE) <<"\n";
            result = false;
        }
        
    }

    return result;
}

OpGraphOp *LISAMapper::getOpNodePtr(OpGraph *opgraph, MRRGNode *n)
{
    int size_vector = opgraph->op_nodes.size();

    OpGraphOp *return_ptr;

    //loop through all op nodes to find same MRRG node as 'this'
    //may be very slow
    for (int i = 0; i < size_vector; i++)
    {
        MRRGNode *current_node = getMappedMRRGNode(opgraph->op_nodes[i]);
        if (current_node == n)
        {
            return_ptr = opgraph->op_nodes[i];
#ifdef DEBUG
            // TODO: check only one occurance of op ptr
#else
            break;
#endif
        }
    }

    return return_ptr;
}

void LISAMapper::mapMRRGNode(OpGraphNode *opnode, MRRGNode *n)
{
    auto mapping_nodes=  mapping[opnode];
    // if(std::find(mapping_nodes.begin(),mapping_nodes.end(), n) !=mapping_nodes.end()){
    //     // std::cout<<"*************find node that has been used\n";
    //     return;
    // }
    try
    {
        mapping[opnode].push_back(n);
    }
    catch (const std::exception &e)
    {
        throw cgrame_error(std::string("LISAMapper Exception Thrown by: [") + e.what() + "] at File: " + std::string(__FILE__) + " Line: " + std::to_string(__LINE__));
    }
    occupancy[n]++;
    // occupancy_detail[n].push_back(opnode);
}

void LISAMapper::mapAllMRRGNodes(OpGraphNode *opnode, std::vector<MRRGNode *> nodes)
{
    // map all nodes
    for (auto &n : nodes)
    {
        mapMRRGNode(opnode, n);
    }
}

void LISAMapper::unmapMRRGNode(OpGraphNode *opnode, MRRGNode *n)
{
    /*
    for(auto & node: mapping[opnode])
    {
        if(node == n)
        {
            //occupancy[n]--;
            //assert(occupancy[n] >= 0);
            mapping[opnode].erase(n);
            break;
        }
    }
    */

    try
    {
        auto iter = find(mapping[opnode].begin(), mapping[opnode].end(), n);
        if (iter != mapping[opnode].end())
        {
            // for(auto it = occupancy_detail[n].begin();it != occupancy_detail[n].end();it++){
            //     if((*it) == opnode){
            //         occupancy_detail[n].erase(it);
            //         break;
            //     }
            // }
            occupancy[n]--;
            assert(occupancy[n] >= 0);
            mapping[opnode].erase(iter);
        }
    }
    catch (const std::exception &e)
    {
        throw cgrame_error(std::string("LISAMapper Exception Thrown by: [") + e.what() + "] at File: " + std::string(__FILE__) + " Line: " + std::to_string(__LINE__));
    }
}

std::vector<MRRGNode *> LISAMapper::unmapAllMRRGNodes(OpGraphNode *opnode)
{
    std::vector<MRRGNode *> result;

    // save mapping
    result = mapping[opnode];

    // unmap all nodes
    for (auto &n : mapping[opnode])
    {
        occupancy[n]--;
        //   for(auto it = occupancy_detail[n].begin();it != occupancy_detail[n].end();){
        //         if((*it) == opnode){
        //             it = occupancy_detail[n].erase(it);
        //             break;
        //         }else{
        //             it ++;
        //         }
        //     }
        assert(occupancy[n] >= 0);
    }

    mapping[opnode].clear();

    return result;
}

/**
  Returns the MRRGNode that the Op is mapped to, NULL if unmapped
 **/
MRRGNode *LISAMapper::getMappedMRRGNode(OpGraphOp *op)
{
    std::vector<MRRGNode *> &mapped_nodes = mapping[op];
    if (mapped_nodes.size() != 0)
    {
        assert(mapped_nodes.size() == 1); // Op should only be mapped to a single node
        assert(mapped_nodes[0]);          // make sure we don't have a NULL pointer

        return mapped_nodes[0];
    }

    return NULL; // return null if unmapped
}

// Route a val node on to the MRRG, false if unable to route

//zhaoying comment: CGRA-ME did not check latency. So the mapping might not be accurate.
// For LISA, we need to map some critical edge first, which needs an accurate latency. It seems we cannot solve this.
// It seems it is not fair to LISA!!!
bool LISAMapper::routeVal(OpGraphVal *val, OpGraphOp *output_op )
{
    assert(val);

    updateOveruse(mrrg_);
    int old_overuse = overuse_num_;
    // cout << "Routing val: " << *val << endl;


#ifdef DEBUG_ROUTING
    cout << "Routing val: " << *val << endl;
    // cout << "Routing from: " << *(val->input->mapped_nodes[0]) << endl;
#endif
    // verify that fanin and fanouts are placed
    if (mapping[val->input].size() == 0)
    {
#ifdef DEBUG_ROUTING
        cout << "Routing input not mappped!" << endl;
#endif
        return false;
    }
    // verify fanouts placed; initialize fanout_mapped flags
    std::map<MRRGNode *, bool> fanout_mapped;
    std::map<MRRGNode *, unsigned int> output_number;

    // TODO: is this the right place to resize this vector?
    val->output_latency.resize(val->output.size());

    for (unsigned int i = 0; i < val->output.size(); ++i)
    {
        OpGraphOp *fos = val->output[i];

        if(output_op != NULL){
            if (fos !=  output_op){
                continue;
            }
        }

        if (mapping[fos].size() == 0)
        {
#ifdef DEBUG_ROUTING
            cout << *fos << "NOT MAPPED?!?!" << endl;
#endif
            assert(0);
            return false;
        }
        // TODO: this next line of code is so ugly.....
        MRRGNode *dst = mapping[fos][0]->operand[val->output_operand[i]];
        #ifdef DEBUG_ROUTING
            std::cout << "dest" << dst->getFullName() << "\n";
        #endif
        output_number[dst] = i;
        fanout_mapped[dst] = false;
        // #ifdef DEBUG_ROUTING
        //         cout << "          to: " << *(fos->mapped_nodes[0]) << ", operand = " << val->output_operand[i] << endl;
        // #endif
    }

    assert(fanout_mapped.size() != 0 );

    std::list<MRRGNode *> src_nodes;
    src_nodes.push_back(mapping[val->input][0]);

    std::map<MRRGNode *, int> min_distance; // if node exites in the map, this node has been mapped.
    min_distance[mapping[val->input][0]] = 0;

//lamda function
    auto find_least_occupy_node = [&](MRRGNode *output_node) {
        int curr_dis = min_distance[output_node] - output_node->latency * 100 - 1;
        int least_occupy_value = occupancy[output_node->prev];
        auto least_occupy_node = output_node->prev;
        // std::cout<<"fan in size "<<output_node->fanin.size()<<"\n";
        for (auto dddd : output_node->fanin)
        {
            if (min_distance.find(dddd) != min_distance.end())
            {
                if (min_distance[dddd] == curr_dis)
                {
                    if (occupancy[dddd] < least_occupy_value)
                    {
                        least_occupy_value = occupancy[dddd];
                        least_occupy_node = dddd;
                    }
                }
            }
        }
        return least_occupy_node;
    };
//lamda function2    
    auto recursive_update = [&](MRRGNode  * node){
        std::queue<MRRGNode *> update_dist;
        update_dist.push(node);
        while(!update_dist.empty()){
            MRRGNode *update_node = update_dist.front();
            update_dist.pop();
            int d_of_updatenode = min_distance[update_node];
            for (auto output_node : update_node->fanout)
            {
                if (min_distance.find(output_node) != min_distance.end())
                {
                    bool need_update = false;
                    if (output_node->prev == update_node)
                    {
                        need_update = true;
                    }
                    else if ((min_distance[output_node] - output_node->latency * 100 ) > d_of_updatenode - 1)
                    {
                        need_update = true;
                    }
                    if (!need_update)
                    {
                        auto least_node = find_least_occupy_node(output_node);
                        output_node->prev = least_node;
                        continue;
                    }

                    if (min_distance[output_node] <= d_of_updatenode + output_node->latency * 100 + 1)
                        continue;
                    output_node->prev = update_node;
                    update_dist.push(output_node);
                    min_distance[output_node] = d_of_updatenode + output_node->latency * 100  + 1;
                }
            }
                    
        }
    };
//lamda function3
    auto update_fanout = [&](MRRGNode *node){
        int min_d = min_distance[node] - node->latency * 100  - 1;
        MRRGNode *min_node = node->prev;
        for (auto node_input : node->fanin)
        {
            if (min_distance.find(node_input) != min_distance.end())
            {
                int d = min_distance[node_input];
                if (d < min_d)
                {
                    min_d = d;
                    min_node = node_input;
                }
            }
        }
        if (min_d == min_distance[node] - node->latency * 100  - 1)
        {
            auto least_node = find_least_occupy_node(node);
            node->prev = least_node;
            return true;
        }
        assert(min_d < min_distance[node] - node->latency * 100  - 1);
// we need to update the value
#ifdef DEBUG_ROUTING
        cout << "-----------update: " << node->getFullName() << " from " << min_distance[node] << " to " << min_d + node->latency << endl;
#endif
        min_distance[node] = min_d + node->latency * 100  + 1;
        assert(min_node != NULL);
        assert(min_node);
        node->prev = min_node;

        recursive_update(node);

        return true;
    };


    bool all_fanouts_mapped = false;
    while (!all_fanouts_mapped)
    {
        // try mapping
        std::priority_queue<std::pair<float, MRRGNode *>, std::vector<std::pair<float, MRRGNode *>>> queue;
        std::vector<MRRGNode *> nodes_in_queue;

        for (auto s = src_nodes.begin(); s != src_nodes.end(); ++s)
        {
#ifdef DEBUG_ROUTING
            cout << "queueing fanouts of src_node: " << **s << endl;
#endif
            for (auto n = (*s)->fanout.begin(); n != (*s)->fanout.end(); ++n)
            {
                // check that fanouts arent already in the src_node list
                // if(occupancy[*n]> 1) continue;

                if (find(src_nodes.begin(), src_nodes.end(), *n) == src_nodes.end() && !(*n)->prev)
                {
                    if (occupancy[*n] >= 1 && ( !allow_overuse || overuse_num_ >=  MAX_OVERUSE ) ){

                        continue;
                    }
                        
                    (*n)->prev = (*s);
                    // std::cout << " add1 fanout" << (*n)->getFullName() << "\n";
                    min_distance[*n] = min_distance[(*s)] + (*n)->latency * 100  + 1;

                    queue.push(std::make_pair(getCost(*n), (*n)));
                    //nodes_in_queue.push_back(*n);
                }
                else
                {
                    // check whether need to update distance
                    // the idea here: when find a  new node which can reach the old node, check the distance of old node and see if need update.
                    // if need update, then we recursively update the fan-in. Meanwhile, for each updated node, let us try to find the best fan-in with the
                    // same least latency but with least occupancy
                    if (!(*n)->prev)
                        continue;
                    if (min_distance.find(*n) == min_distance.end())
                        continue;
                    update_fanout(*n);
                }
            }
        }

        // while we still have more routing options, try to route
        bool mapped_a_node = false;
        while (!queue.empty())
        {
            auto node = queue.top();
            queue.pop();
            //nodes_in_queue.erase(find(nodes_in_queue.begin(), nodes_in_queue.end(), node.second));

            float node_cost = node.first;
            MRRGNode *n = node.second;
#ifdef DEBUG_ROUTING
            cout << "queue popped: " << *n << ": " << node_cost << endl;
#endif
            // if n is a new sink
            // test if n->first is in map, then test if is unmapped
            if (fanout_mapped.find(n) != fanout_mapped.end() && !fanout_mapped[n])
            {
#ifdef DEBUG_ROUTING
                cout << "Router found sink: " << *n << "." << endl;
#endif
                // backtrack and mark path
                // FUTURE TODO: mark dst OpGraphOp node with path statistics: delay, cycle latency etc..

                // add dst node (the node that is the FU input) to the src list.
                src_nodes.push_back(n);

                unsigned int latency = 0;
                MRRGNode *k = n->prev;

                while (find(src_nodes.begin(), src_nodes.end(), k) == src_nodes.end())
                {
#ifdef DEBUG_ROUTING
                    cout << *k << endl;
#endif
                    //if((*k)->type == MRRG_NODE_REGISTER)
                    //{
                    latency += k->latency;
                    //}

                    // add node to src node list
                    src_nodes.push_back(k);
                    // backtrack
                    k = k->prev;
                }

                // clear queue and unmark all backtrack paths to src_nodes
                // inner while loop will end when queue is empty
                while (!queue.empty())
                {
                    auto node_pair = queue.top();
                    queue.pop();
                    //nodes_in_queue.erase(find(nodes_in_queue.begin(), nodes_in_queue.end(), node_pair.second));

                    MRRGNode *k = node_pair.second;
                    // while we haven't backtracked to already removed path AND we havent reached any source nodes
                    while (k != NULL && find(src_nodes.begin(), src_nodes.end(), k) == src_nodes.end())
                    {
                        MRRGNode *prev = k;
                        k = k->prev;
                        prev->prev = NULL;
                        auto it = min_distance.find(prev);
                        if (it != min_distance.end())
                            min_distance.erase(it);
                    }
                }
                for (auto it = min_distance.begin(); it != min_distance.end();)
                {
                    auto node = it->first;
                    if (find(src_nodes.begin(), src_nodes.end(), node) == src_nodes.end())
                    {
                        it = min_distance.erase(it);
                    }
                    else
                    {
                        it++;
                    }
                }

                // we have now mapped the fanout
                // std::cout<<"mapped one fanout\n";
                #ifdef DEBUG_ROUTING
                        cout << "latency "<<n->getFullName()<<" "<<latency << endl;
                #endif
                fanout_mapped[n] = true;
                val->output_latency[output_number[n]] = latency;
                mapped_a_node = true;
            }
            else
            {
                // queue up all fanouts
                if (n->type == MRRG_NODE_ROUTING)
                {
                    for (auto i = n->fanout.begin(); i != n->fanout.end(); ++i)
                    {
                        // check that fanouts arent already in the src_node
                        //  if(occupancy[*i]> 1) continue;
                        if (!(*i)->prev)
                        {
                            if (occupancy[*i] >= 1 && ( !allow_overuse || overuse_num_ >=  MAX_OVERUSE ) ){
                                continue;
                            }
                            (*i)->prev = n;
                            // std::cout << " add2 fanout" << (*i)->getFullName() << "\n";

                            queue.push(std::make_pair(node_cost + getCost(*i), (*i)));
                            min_distance[*i] = min_distance[n] + (*i)->latency  * 100 + 1;
                        }
                        else
                        {
                            // check whether need to update distance
                            if (min_distance.find(*i) == min_distance.end())
                                continue;
                            update_fanout(*i);
                        }
                    }
                }
            }
        }

        // if we exited the loop and didn't map a node, we have an unreachable route
        if (!mapped_a_node) // if we exited the loop and didn't map a node, we have an unreachable route
        {
#ifdef DEBUG_ROUTING
            cout << "Routing Failed" << endl;
#endif
            return false;
        }
        all_fanouts_mapped = true;
        for (auto i = fanout_mapped.begin(); i != fanout_mapped.end(); ++i)
        {
            all_fanouts_mapped &= i->second;
        }
       
    }
   

    for (auto s = src_nodes.begin(); s != src_nodes.end(); ++s)
    {
        if ((*s)->type == MRRG_NODE_ROUTING)
        {
            #ifdef DEBUG_ROUTING
            cout << "Mapping routing node: " << **s << "  "<<(*s)->latency<<"\n";
            #endif
            // map val to this node
            mapMRRGNode(val, *s);
        }
        else
        {
#ifdef DEBUG_ROUTING
            cout << "non routing node in source node list" << endl;
#endif
        }
    }
  
#ifdef DEBUG_ROUTING
    cout << "Routing Succeeded" << endl;
#endif

   
    min_distance.clear();

     updateOveruse(mrrg_);
    int current_overuse = overuse_num_;
    int overuse_diff = current_overuse - old_overuse;
    overuse_counter[val->input] =overuse_counter[val->input]  + overuse_diff;
    for(auto output_node: val->output){
        overuse_counter[output_node] =overuse_counter[output_node]  + overuse_diff;
    }
    // assert(false);
    return true;
}

float LISAMapper::getTotalOpCost(OpGraphOp *op)
{
    float result = 0.0;
    result += getCost(op);

    if (op->opcode != OPGRAPH_OP_INPUT && op->opcode != OPGRAPH_OP_CONST)
    {
        for (auto in = op->input.begin(); in != op->input.end(); ++in)
        {
            result += getCost((*in));
        }
    }

    if (op->opcode != OPGRAPH_OP_OUTPUT && op->opcode != OPGRAPH_OP_STORE)
    {
        result += getCost(op->output);
    }

    

    return result;
}

// Rips up Op as well as routes in and out (if they exist) and records the cost of whatever was ripped up
// if a pointer to cost is given, the cost is also returned
OpMapping LISAMapper::ripUpOp(OpGraphOp *op, float *cost)
{
    OpMapping result;

    int id = node_to_id_[op];
    if(mapping[op].size()>0){
        auto old_fu = mapping[op].front();
        mpos mmmm{old_fu->x_, old_fu->y_, old_fu->cycle};
            
        for(auto it =  mapped_pos.begin(); it!= mapped_pos.end(); ){
            if(it->x == mmmm.x && mmmm.y == it->y  && mmmm.t == it->t){
                it = mapped_pos.erase(it);
            }else{
                it ++;
            }
        }
    }

    if (cost)
    {
        *cost = 0.0;
        *cost += getCost(op);
    }

    result.mapping[op] = unmapAllMRRGNodes(op);

    if (op->opcode != OPGRAPH_OP_INPUT && op->opcode != OPGRAPH_OP_CONST)
    {
        for (auto &in : op->input)
        {
            if (cost)
            {
                *cost += getCost(in);
            }
            // record input val mapping, only if not already recorded. This happens when the same value feeds both operands.
            if (result.mapping.find(in) == result.mapping.end())
                result.mapping[in] = unmapAllMRRGNodes(in);
        }
    }

    if (op->opcode != OPGRAPH_OP_OUTPUT && op->opcode != OPGRAPH_OP_STORE)
    {
        if (cost)
        {
            *cost += getCost(op->output);
        }

        // record output val mapping
        result.mapping[op->output] = unmapAllMRRGNodes(op->output);
    }

    return result;
}

void LISAMapper::restoreOp(OpMapping oldmap)
{
    for (auto m = oldmap.mapping.begin(); m != oldmap.mapping.end(); ++m)
    {
        mapAllMRRGNodes(m->first, m->second);
    }
}

// get a random unoccupied  FU
MRRGNode *LISAMapper::getRandomUnoccupiedFU(MRRG *mrrg, OpGraphOp *op)
{
    std::vector<MRRGNode *> candidates;
    for (auto &fu : mrrg->function_nodes)
    {
        if (fu->canMapOp(op) && occupancy[fu] == 0)
        {
            candidates.push_back(fu);
        }
    }

    if (candidates.size() < 1)
    {
        cout << "Could not place: " << *op << endl;
        assert(candidates.size() > 0);
    }

    return candidates[rand() % candidates.size()];
}
/*
bool compare_mrrg_node_occupancy(MRRGNode* a, MRRGNode* b)
{
    return a->occupancy < b->occupancy;
}

MRRGNode* getRandLeastOccupied(std::vector<MRRGNode*> nodes)
{
    assert(nodes.size() > 0);
    std::random_shuffle(nodes.begin(), nodes.end());
    std::sort(nodes.begin(), nodes.end(), compare_mrrg_node_occupancy);

    return nodes[0];
}
*/

/*
bool compare_pair_mrrg_node_occupancy(std::pair<MRRGNode*,int> a, std::pair<MRRGNode*,int> b)
{
    return a.first->occupancy < b.first->occupancy;
}

bool compare_pair_mrrg_node_dist(std::pair<MRRGNode*,int> a, std::pair<MRRGNode*,int> b)
{
    return a.second < b.second;
}

MRRGNode* getRandClosestLeastOccupied(std::vector<std::pair<MRRGNode*, int> > nodes)
{
    assert(nodes.size() > 0);
    std::random_shuffle(nodes.begin(), nodes.end());
    std::sort(nodes.begin(), nodes.end(), compare_pair_mrrg_node_occupancy);

    std::vector<std::pair<MRRGNode*, int>>::iterator firstgroupend;

    int initial_occupancy = nodes[0].first->occupancy;
    std::vector<std::pair<MRRGNode*, int> > nodes2;

    for(auto fu = nodes.begin(); fu != nodes.end(); fu++)
    {
        if((*fu).first->occupancy != initial_occupancy)
        {
            firstgroupend = fu;
            break;
        }
        else
        {
            nodes2.push_back(*fu);
        }

    }

    std::sort(nodes2.begin(), nodes2.end(), compare_pair_mrrg_node_dist);

    return nodes2[rand() % nodes2.size()].first;
}
*/

/*
MRRGNode* getInputMRRGNode(OpGraphOp* op, int operand)
{
    if(op->opcode == OPGRAPH_OP_INPUT)
    {
        return NULL;
    }

    if(operand < op->input.size())
    {
        OpGraphVal* input_val = op->input[operand];

        assert(input_val->input);
        return getMappedMRRGNode(input_val->input);
    }
    else
    {
        return NULL;
    }
}
MRRGNode* getNearbyFU(MRRG* mrrg, OpGraphOp* op)
{
    // Check input ops to see if they are placed
    MRRGNode* f0 = getInputMRRGNode(op, 0);
    MRRGNode* f1 = getInputMRRGNode(op, 1);

    // intersect f0's nearby fu's and f1's nearby fu's
    std::vector<std::pair<MRRGNode*, int> > candidates;
    if(f0 && f1)
    {
        candidates = candidate_fu_intersect(f0->neighbourFUs, f1->neighbourFUs);
    }
    else if(f0)
    {
        candidates = f0->neighbourFUs;
    }
    else if(f1)
    {
        candidates = f1->neighbourFUs;
    }
    // else empty set...
    // maybe try and match output?
    else
    {

    }

    // and filter FU's that cant map the OP, then choose the closest and least congested
    std::vector<std::pair<MRRGNode*, int> > filtered_candidates;
    for(auto & fu : candidates)
    {
        if(fu.first->canMapOp(op))
        {
            filtered_candidates.push_back(fu);
        }
    }

    if(filtered_candidates.size() < 1)
    {
#ifdef PLACEMENT_DEBUG
        cout << "Could not getNearbyFU for : " << *op << ". Revert to getRandomFU()" << endl;
#endif
        // assert(filtered_candidates.size() > 0);
        return getRandomFU(mrrg, op);
    }

    return getRandClosestLeastOccupied(filtered_candidates);
}
*/
bool LISAMapper::placeOp(OpGraphOp *op, MRRGNode *n)
{
    assert(n->type == MRRG_NODE_FUNCTION);

    if (mapping[op].size() != 0) // Make sure this op is not placed anywhere else
        return false;

    mapMRRGNode(op, n);

    mapped_pos.insert(mpos{n->x_, n->y_, n->cycle});

    return true;
}

inline bool accept(float delta_cost, float temperature)
{
    if (delta_cost < 0)
        return true;

    float probability = exp(-(delta_cost) / temperature);

    return probability > ((float)rand() / (float)RAND_MAX);
}

inline bool accept(float new_cost, float old_cost, float temperature)
{
    if (new_cost < old_cost)
        return true;

    float probability = exp(-(new_cost - old_cost) / temperature);

    return probability > ((float)rand() / (float)RAND_MAX);
}

inline float LISAMapper::updateTempConst(float t)
{
    return t * const_temp_factor; //0.999;
}

bool LISAMapper::routeOp_withfailed_node(OpGraphOp *op, MRRG *mrrg, std::set<OpGraphOp *> & failed)
{
    bool result = true;
    // route input vals
    if (op->opcode != OPGRAPH_OP_INPUT && op->opcode != OPGRAPH_OP_CONST)
    {
        for (auto v = op->input.begin(); v != op->input.end(); ++v)
        {
            // check all nodes and make sure prev is NULL
            for (auto node = mrrg->function_nodes.begin(); node != mrrg->function_nodes.end(); ++node)
                (*node)->prev = NULL;

            for (auto node = mrrg->routing_nodes.begin(); node != mrrg->routing_nodes.end(); ++node)
                (*node)->prev = NULL;

            bool r = routeVal(*v);
            updateOveruse(mrrg);

            if (!r)
            {   
                failed.insert((*v)->input);
                // cout << "Could not route val: \"" << **v << "\". Turn on DEBUG_ROUTING for more info." << endl;
            }
            result &= r;
        }
    }

    // route output val
    if (op->opcode != OPGRAPH_OP_OUTPUT && op->opcode != OPGRAPH_OP_STORE)
    {
        // check all nodes and make sure prev is NULL
        for (auto node = mrrg->function_nodes.begin(); node != mrrg->function_nodes.end(); ++node)
            (*node)->prev = NULL;

        for (auto node = mrrg->routing_nodes.begin(); node != mrrg->routing_nodes.end(); ++node)
            (*node)->prev = NULL;

        bool r = routeVal(op->output);
        updateOveruse(mrrg);
        if (!r)
        {
            for(auto o: op->output->output){
                failed.insert(o);
            }
            
            // cout << "Could not route val: \"" << *(op->output) << "\". Turn on DEBUG_ROUTING for more info." << endl;
        }
        result &= r;
    }

    return result;
}
bool LISAMapper::routeOp(OpGraphOp *op, MRRG *mrrg)
{
    bool result = true;
    // route input vals
    if (op->opcode != OPGRAPH_OP_INPUT && op->opcode != OPGRAPH_OP_CONST)
    {
        for (auto v = op->input.begin(); v != op->input.end(); ++v)
        {

            if(!(*v)) continue;
            // check all nodes and make sure prev is NULL
            for (auto node = mrrg->function_nodes.begin(); node != mrrg->function_nodes.end(); ++node)
                (*node)->prev = NULL;

            for (auto node = mrrg->routing_nodes.begin(); node != mrrg->routing_nodes.end(); ++node)
                (*node)->prev = NULL;

            bool r = routeVal(*v);
            updateOveruse(mrrg);
            if (!r)
            {
                // cout << "Could not route val: \"" << **v << "\". Turn on DEBUG_ROUTING for more info." << endl;
            }
            result &= r;
        }
    }

    // route output val
    if (op->opcode != OPGRAPH_OP_OUTPUT && op->opcode != OPGRAPH_OP_STORE)
    {
        // check all nodes and make sure prev is NULL
        for (auto node = mrrg->function_nodes.begin(); node != mrrg->function_nodes.end(); ++node)
            (*node)->prev = NULL;

        for (auto node = mrrg->routing_nodes.begin(); node != mrrg->routing_nodes.end(); ++node)
            (*node)->prev = NULL;

        bool r = routeVal(op->output);
        updateOveruse(mrrg);
        if (!r)
        {
            // cout << "Could not route val: \"" << *(op->output) << "\". Turn on DEBUG_ROUTING for more info." << endl;
        }
        result &= r;
    }

    return result;
}

bool LISAMapper::routeOpInput(OpGraphOp *op, MRRG *mrrg)
{
    bool result = true;
    // route input vals
    if (op->opcode != OPGRAPH_OP_INPUT && op->opcode != OPGRAPH_OP_CONST)
    {
        for (auto v = op->input.begin(); v != op->input.end(); ++v)
        {

            if(!(*v)) continue;
            // check all nodes and make sure prev is NULL
            for (auto node = mrrg->function_nodes.begin(); node != mrrg->function_nodes.end(); ++node)
                (*node)->prev = NULL;

            for (auto node = mrrg->routing_nodes.begin(); node != mrrg->routing_nodes.end(); ++node)
                (*node)->prev = NULL;

            bool r = routeVal(*v, op);
            updateOveruse(mrrg);
            if (!r)
            {
                // cout << "Could not route val: \"" << **v << "\". Turn on DEBUG_ROUTING for more info." << endl;
            }
            result &= r;
        }
    }

    // route output val


    return result;
}

static inline double getcurrenttime()
{
    struct timeval t;

    gettimeofday(&t, NULL);

    return t.tv_sec + t.tv_usec * 0.000001;
}

// This is the main mapping function
// true on success, false on failure
Mapping LISAMapper::mapOpGraph(std::shared_ptr<OpGraph> opgraph, int II, std::string arch_model_name)

{
     for(int i = 0; i< systolic_array_x; i++){
        for(int j=0; j< systolic_array_y;j++){
            systolic_pe_index.emplace(std::make_pair(i,j),sys_arr[i][j]);
        }
    }
    arch_model_name_ = arch_model_name;

    for(auto op: opgraph->op_nodes){
        overuse_counter.emplace(op, 0);
    }

    // init op_edges
    {
        for(auto node: opgraph->op_nodes){
            if (node->opcode != OPGRAPH_OP_INPUT && node->opcode != OPGRAPH_OP_CONST)
            {   
                int node_id = node_to_id_[node];
                assert(dfg_label_->find(node_id)!= dfg_label_->end());
                auto & ass_value = dfg_label_->at(node_id).association;
                for (auto input_value: node->input){
                    if(input_value == 0) continue;
                    auto input_node = input_value->input;
                    int input_node_id = node_to_id_[input_node];
                    if(node_id ==  input_node_id) continue;
                    assert(ass_value.find(input_node_id) != ass_value.end());
                    std::cout<<"add edge"<<input_node->name <<"->"<<node->name<<"\n";
                    op_edges_.emplace(op_edge{input_node, node, ass_value[input_node_id].second});
                }
            }
        }
    }

    
    mapper_start_time  = std::chrono::steady_clock::now();
    this->II_ = II;
    opgraph_ = opgraph;
    MAX_OVERUSE = opgraph->op_nodes.size() * 100;
    std::cout<<"cgra size ("<<cgra->ROWS<<", "<<cgra->COLS<<")"<<std::endl;
    // set_edges(opgraph);
    // find_backedge();
    std::set<int> node_list;
    for(auto node: id_to_node_){
        node_list.insert(node.first);
    }

    
    
    // assert(false);
    // return NULL;
    // get the mrrg object
    MRRG *mrrg = cgra->getMRRG(II).get();
    mrrg_ = mrrg;
     if (arch_model_name_.find("systolic") != std::string::npos) {
        mrrg_->makeSystolicArray(systolic_pe_index, systolic_array_x);
    }

    
     if (arch_model_name_.find("leftmostmemory") != std::string::npos) {
        mrrg_->makeLeftMostMemoryAccess();
    }
    
    for(auto node: mrrg->routing_nodes){
        occupancy[node] = 0;
    }
    for(auto node: mrrg->function_nodes){
        occupancy[node] = 0;
    }
    // Create result obj
    Mapping mapping_result(cgra, II, opgraph);

    // Set the random seed
    srand(this->rand_seed);

#ifdef ANNEAL_DEBUG
    ofstream anneal_debug;
    anneal_debug.open("anneal_debug.csv");
#endif

    // Sort the opgraph nodes
    // TODO: This does not work on graphs with back edges
    // topological_sort(opgraph);

    cout << "using LISA" << std::endl;

    LOG(LISAGNN) << "Initial placement:" << "\n";
    // randomized initial placement

    std::vector<OpGraphOp*>  sorted_nodes;
    for(auto node: opgraph->op_nodes){
        if(find(sorted_nodes.begin(), sorted_nodes.end(),node) == sorted_nodes.end()){
            sorted_nodes.push_back(node);
        }
    }
    if(!lisa_eval_routing_priority){
        std::sort (sorted_nodes.begin(), sorted_nodes.end(), 
        [&, this](OpGraphOp* a, OpGraphOp* b) {
            int a_order = dfg_label_->at(node_to_id_[a]).schedule_order;
            int b_order = dfg_label_->at(node_to_id_[b]).schedule_order;
            if(a_order!= b_order) return a_order < b_order;
            int a_id = node_to_id_[a];
            int b_id = node_to_id_[b];
            // if(node_children_[a_id].find(b_id) != node_children_[a_id].end()) return true;
            // if(node_children_[b_id].find(a_id) != node_children_[b_id].end()) return false;
            return node_asap_[a_id] < node_asap_[b_id]; 
        });
    }

    std::vector<OpGraphOp*> sorted_nodes_by_edge;
    for(auto edge: op_edges_){
        auto src= edge.src;
        auto des = edge.des;
        std::cout<<"edge: "<<src->name<<" -> "<<des->name<<"\n";

        auto src_iter = std::find(sorted_nodes_by_edge.begin(), sorted_nodes_by_edge.end(), src);
        
        if (src_iter == sorted_nodes_by_edge.end()){
            sorted_nodes_by_edge.push_back(src);
        }
        auto des_iter = std::find(sorted_nodes_by_edge.begin(), sorted_nodes_by_edge.end(), des);
            if (des_iter == sorted_nodes_by_edge.end()){
            sorted_nodes_by_edge.push_back(des);
        }
        // std::cout<<"iter:"<< i<<" src:"<<edge.src->name<<" des:"<<edge.des->name<<" route length "<<edge.length<<"\n";
        // routed = routeVal(edge.src->output, edge.des);
        //  if(!routed) break;
    }
    std::cout<<"sorted nodes: ";
    for(auto op: sorted_nodes){
        std::cout<< op->name<<"  ";
    }
    std::cout<<"\ntemp_sored_nodes : ";
    for(auto op: sorted_nodes_by_edge){
        std::cout<< op->name<<"  ";
    }
    assert(sorted_nodes.size() == sorted_nodes_by_edge.size());


    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_int_distribution<int> uniform_dist(0, sorted_nodes.size()-1);
       
    bool  routed = true; 
    std::cout<<"sorted nodes\n";
     for (auto &op : sorted_nodes)
        {
            std::cout<< *op<<" id:"<<node_to_id_[op]<< std::endl;
        }
    std::map<OpGraphNode*, MRRGNode*>  bestPlacement;
    int min_cost  = 10000000;
    for(int i  = 0 ; i < 100 ; i++){
        routed = true;
        
         std::map<OpGraphNode*, MRRGNode*>  currPlacement;
         if(place_and_routing_parallel){
            for (auto &op : sorted_nodes)
            {
                MRRGNode *fu;
                if(!lisa_eval_routing_priority){ 
                     fu = getLISAFU(opgraph, mrrg, op);
                }else{
                     fu = getRandomFU(mrrg, op);

                }
                LOG(LISAFU)<<"iter:"<< i <<" " << *op << " : " << *fu ;
                bool placed = placeOp(op, fu);
                currPlacement.emplace(op, fu);
                
                assert(placed);
                routed = routeOpInput(op, mrrg);
                if(!routed) break;
            }
         }else{
            for (auto &op : sorted_nodes)
                {
                    MRRGNode *fu;

                    if(!lisa_eval_routing_priority){ 
                        fu = getLISAFU(opgraph, mrrg, op);
                    }else{
                        fu = getRandomFU(mrrg, op);

                    }
                    LOG(LISAFU)<<"iter:"<< i <<" " << *op << " : " << *fu <<"\n";
                    bool placed = placeOp(op, fu);
                    currPlacement.emplace(op, fu);
                    
                    assert(placed);
                }

                // for (auto &op :sorted_nodes)
                // {
                //     routed = routeOp(op, mrrg);
                    
                //     if(!routed) break;
                // }
               
                 for (auto &op :sorted_nodes_by_edge)
                {
                    routed = routeOp(op, mrrg);
                    
                    if(!routed) break;
                }

         }
        

        

        if(routed){
            if(is_training){
                break;
            }else{
                int curr_cost = getCost(mrrg);
                if(curr_cost < min_cost){
                    bestPlacement =  currPlacement;
                    min_cost = curr_cost;
                }
            }
        } 
        for (auto &op : sorted_nodes)
        {
            ripUpOp(op);
        }
        

        if( i == 50) {
            allow_overuse =  true;
            LOG(LISAGNN) <<"iter:"<< i <<" set overuse\n" ;
        }
        
    }
    if(!routed) return mapping_result;
    if(!is_training){
         for (auto &op : sorted_nodes)
        {
            ripUpOp(op);
        }
          for (auto &op : sorted_nodes)
        {
            placeOp(op, bestPlacement[op]);
        }
        for (auto &op :sorted_nodes_by_edge)
        {
            routed = routeOp(op, mrrg);
            
            assert(routed);
        }
    }
    

    allow_overuse =  false;
    finish_init = true;

    
    // Verify initial place and route
    for (auto &op : opgraph->op_nodes)
    {
        if (mapping[op].size() == 0)
        {
            cout << "No initial placement for: " << op->name << endl;
            assert(0);
        }
    }

    for (auto &val : opgraph->val_nodes)
    {
        if (mapping[val].size() == 0)
        {
            cout << "No initial routing for: " << val->name << endl;
            assert(0);
        }
    }
  

    for (auto &node : mrrg->routing_nodes)
    {
        if (occupancy[node] > node->capacity)
        {
            
            LOG(OVERUSE) << *node << " is overused. (" << occupancy[node] << "/" << node->capacity << ")\n";
            // for(auto mappinged_node: occupancy_detail[node]){
            //     LOG(OVERUSE) <<"use "<<mappinged_node->name;
            // }
            //  LOG(OVERUSE) <<"\n";
        }
        
    }
    // assert(false);
    
    LOG(DUMPMAPPING)<<printMapping()<<"\n";
    // assert(false);
    

#ifdef CALCULATE_INITIAL_TEMPERATURE
    /************************* TRY TO DO 100 iterations before annealing********/

    cout << "Finding delta Costs:" << endl
         << "INITIAL:" << endl
         << "OpGraph Cost: " << getCost(opgraph.get()) << endl;
    cout << "MRRG Cost: " << getCost(mrrg) << endl;
    cout << "Penalty Factor: " << pfactor << endl;

    float max_delta_cost = 0;
    float old_cost = getCost(mrrg);

    //  assert(false);

    //do 100x of perturbation
    int max_try = 100;
    for (int i = 0; i < 100; i++) 
    // for (int i = 0; i < 100; i++)
    {
        //first get a random index
        int vector_size = opgraph->op_nodes.size();
        int index = rand() % vector_size;

        //perturb at this index
        OpGraphOp *op = opgraph->op_nodes[index];

        //substitute a new fu
        MRRGNode * old_fu = mapping[op].front();
        MRRGNode *fu;
        fu = getRandomFU(mrrg, op);
    
        

        if (fu->canMapOp(op))
        {
            //check if occupied
            if (occupancy[fu] == 0)
            {
                //move there
                OpMapping oldmap = ripUpOp(op);
                bool success = placeOp(op, fu);
                if (success)
                    success = routeOp(op, mrrg);
                if(!success ){
                    ripUpOp(op);
                    restoreOp(oldmap);
                    continue;
                } 

                //get new cost
                float new_cost = getCost(mrrg);
                float change_cost = abs(new_cost - old_cost);
                if (change_cost > max_delta_cost)
                    max_delta_cost = change_cost;

                //restore changes
                ripUpOp(op);
                restoreOp(oldmap);
            }
            else
            {
                assert(false);
            }
        }
        else
        {
            //cant map, skip
            i--;
        }
    }

    //calculate initial temperature
    float accept_percentage = 0.99;
    float natural_log = log(accept_percentage);

    // initial temperature
    float temperature = (-1) * max_delta_cost / (natural_log);
    cout << " max delta cost is " << max_delta_cost << endl;
#else
    float temperature = 1000000.0;
#endif
    cout << "Initial Temperature is " << temperature << endl;

    /************************** Done 100 iterations **************************************/

    cout << "Begin annealing" << endl;
    bool no_timelimit = (timelimit == 0.0);
    double start_time = getcurrenttime();
    double current_time = getcurrenttime();

    float current_cost = getCost(mrrg);
    while (no_timelimit || (current_time - start_time) < timelimit)
    {   
        
        // dfg_label = lisa_ctl->getCurrLabel();
        float accept_rate = 0.0;
        float previous_cost = current_cost;
        LOG(LISAGNN) << "########################################\n Annealing at:"
                << "\ttemp: " << temperature
                << "\tpfactor: " << pfactor
                << "\tcost:"<<current_cost ;
        if (inner_place_and_route_loop(opgraph.get(), mrrg, temperature, &accept_rate))
        {
            mapping_result.setMapping(mapping);
            LOG(LISAGNN) << "MappingTime: " << (int)(getcurrenttime() - start_time);
            LOG(LISAGNN) << "MapperTimeout: 0";
            LOG(LISAGNN) << "Mapped: 1";
            mapping_result.setMapped(true);
              if(is_training){
                    current_cost = getCost(mrrg);
                    optimizeMapping(opgraph.get(), mrrg);
                    std::cout<<"lisa training optimize mapping cost from "<<current_cost << " to "<<getCost(mrrg)<<"\n";
                }
            return mapping_result;
        }
        current_cost = getCost(mrrg);

        LOG(LISAGNN) << "aftter SA -> temp:"<<temperature << ","
        << "cost:"<<current_cost << ","
        << pfactor << ","
        << "\n";
        // TODO: Changed from 0.01
        //if(temperature  < 0.001 * mrrg->getCost(pfactor) / mrrg->routing_nodes.size())
        if (accept_rate < cold_accept_rate && current_cost >= previous_cost)
        {
#ifdef ANNEAL_DEBUG
            anneal_debug.close();
#endif
            LOG(LISAGNN) << "Mapper is Cold and no valid mapping was found."<<
            "Current Temperature acceptance rate was: " << accept_rate
            << "Cold acceptance rate is: " << cold_accept_rate
            << "Current Cost is: " << current_cost
            << "Previous  Cost was: " << previous_cost
            << "Mapped: 0";
            return mapping_result;
        }
        LOG(LISAGNN) << "mrrg cost: " << getCost(mrrg)
        << " mrrg size: " << mrrg->routing_nodes.size()
        << " update temp. & pfactor";
        // update temperature
#ifdef SIMPLE_ANNEAL_SCHEDULE
        temperature = updateTempConst(temperature);
#else
        temperature = updateTemperature(temperature, (float)total_accepted / total_tries);
#endif
        // update overuse penalty
        pfactor = pfactor * pfactor_factor;

        current_time = getcurrenttime();

        LOG(LISAGNN) << "current run time: " << (int)(current_time - start_time);

    }
  
    cout << "MapperTimeout: 1" << endl;
    cout << "Mapped: 0" << endl;

    return mapping_result;
}
bool LISAMapper::optimizeMapping(OpGraph *opgraph, MRRG *mrrg){
    int num_swaps = 1000;
    int total_accepted = 0;
    for (int i = 0; i < num_swaps; i++)
    {
        // Get an op
        OpGraphOp *op = opgraph->op_nodes[rand() % opgraph->op_nodes.size()];

        // Get an fu
        MRRGNode *old_fu = mapping[op].front();
        // MRRGNode *fu = getRandomFU(mrrg, op);
        //shared_ptr<OpGraph > my_ptr(opgraph);
        MRRGNode *fu;
        fu = getCloseRandomFU( mrrg, op, old_fu, 2, 2);
           
         
        // std::cout<<" old fu "<<old_fu->getFullName()<<" new fu: "<<fu->getFullName()<<"\n";
        // make sure that it is a different FU
        if (fu == getMappedMRRGNode(op))
            continue;

        //find the cost of this selected op
        float old_cost = getCost(mrrg); //getTotalOpCost(op);
        
        if (fu->canMapOp(op))
        {
            //check if occupied
            if (occupancy[fu] == 0)
            {
                //move there
                OpMapping oldmap = ripUpOp(op);
                bool success = placeOp(op, fu);
                if (!success)
                {
                    cout << "Could not place OP" << endl;
                }

                if (success)
                {
                    success = routeOp(op, mrrg);
                    if (!success)
                    {
                        // cout << "Could not route OP. It is likely that there is a disconnect in the architecture." << endl;
                        // assert(success);
                        ripUpOp(op);
                        restoreOp(oldmap);
                        continue;
                    }
                }

                //get new cost
                float new_cost_x = getCost(mrrg); // getTotalOpCost(op);
                float delta_cost = new_cost_x - old_cost;

                if (accept(delta_cost, 0))
                {
                    total_accepted++;
                }
                else
                {
                    //restore changes
                    ripUpOp(op);
                    restoreOp(oldmap);
                }
            }
            else // there must only be one unit mapped here
            {
                assert(false);
            }
        }
    }
     bool mrrg_overuse = checkOveruse(mrrg);
     assert(mrrg_overuse);
    return true;

    

}

bool LISAMapper::inner_place_and_route_loop(OpGraph *opgraph, MRRG *mrrg, float temperature, float *accept_rate)
{
    std::cout<<"swap_factor:"<<swap_factor<<"\n";
    int num_swaps = swap_factor ;
    if(is_training) num_swaps = swap_factor  ;
    int total_accepted = 0;
    int total_tries = 0;

    for (int i = 0; i < num_swaps; i++)
    {
        if(is_training){
             auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-mapper_start_time;
            int used_second =   elapsed_seconds.count() ;
            if(used_second > 3600* 3) return false;
        }
       
        // Get an op
        int node_max_overuse = 0, node_min_overuse = 99999999;
        // std::cout<<"overuse counter:";
        for(auto overuse_pair: overuse_counter){
            node_max_overuse = std::max(node_max_overuse, overuse_pair.second);
            node_min_overuse = std::min(node_min_overuse, overuse_pair.second);
            // std::cout<<overuse_pair.first->name<<", "<<overuse_pair.second<<"  ";
        }
        int overuse_diff = node_max_overuse - node_min_overuse;
        // std::cout<<"\n";

        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> d{0,3};
        int selected_id = std::abs(std::round(d(gen)));
        auto op_sort_by_overuse = opgraph->op_nodes;
        std::sort(op_sort_by_overuse.begin(), op_sort_by_overuse.end(),[&](OpGraphOp * a, OpGraphOp * b){
                return overuse_counter[a] > overuse_counter[ b];});
        // if(selected_cost > max_diff)  selected_cost = max_diff;
        // selected_cost  = min_cost + selected_cost;
        // Generator g(0, 1000, -overuse_diff, overuse_diff);
        // std::cout<<"selected overuse number "<<selected_id<<"\n";
        if(selected_id >=  op_sort_by_overuse.size()) selected_id =   op_sort_by_overuse.size() - 1;
        // std::cout<<"selected op"<<op_sort_by_overuse[selected_id]->name<<"\n";


        // OpGraphOp *op = opgraph->op_nodes[rand() % opgraph->op_nodes.size()];
        OpGraphOp *op ;
        // if(is_training){
            op = opgraph->op_nodes[rand() % opgraph->op_nodes.size()];
        // }else{
            // op = op_sort_by_overuse[selected_id];
        // }
        // Get an fu
        MRRGNode *old_fu = mapping[op].front();
        // MRRGNode *fu = getRandomFU(mrrg, op);
        //shared_ptr<OpGraph > my_ptr(opgraph);
        MRRGNode *fu;
        if (is_training){
            fu = getCloseRandomFU( mrrg, op, old_fu);
           
        }else  if(lisa_eval_routing_priority){ 
            fu = getRandomFU(mrrg, op);
        }else{
            // fu = getCloseRandomFU( mrrg, op, old_fu,1 ,1);
           fu = getLISAFU(opgraph_, mrrg, op, total_accepted, total_tries, num_swaps);
        }
         
        // std::cout<<" old fu "<<old_fu->getFullName()<<" new fu: "<<fu->getFullName()<<"\n";
        // make sure that it is a different FU
        if (fu == getMappedMRRGNode(op))
            continue;

        //find the cost of this selected op
        float old_cost = getCost(mrrg); //getTotalOpCost(op);
        
        if (fu->canMapOp(op))
        {
            // This is an actual attempt to swap
            total_tries++;

            //check if occupied
            if (occupancy[fu] == 0)
            {
                //move there
                OpMapping oldmap = ripUpOp(op);
                bool success = placeOp(op, fu);
                if (!success)
                {
                    cout << "Could not place OP" << endl;
                }

                if (success)
                {
                    success = routeOp(op, mrrg);
                    if (!success)
                    {
                        // cout << "Could not route OP. It is likely that there is a disconnect in the architecture." << endl;
                        // assert(success);
                        ripUpOp(op);
                        restoreOp(oldmap);
                        continue;
                    }
                }

                //get new cost
                float new_cost_x = getCost(mrrg); // getTotalOpCost(op);
                float delta_cost = new_cost_x - old_cost;

                if (accept(delta_cost, temperature))
                {
                    total_accepted++;
                }
                else
                {
                    //restore changes
                    ripUpOp(op);
                    restoreOp(oldmap);
                }
            }
            else // there must only be one unit mapped here
            {
#ifndef ALLOW_MULTIPLE_PLACEMENT
                assert(fu->occupancy == 1);
#endif

                //swap, rip off two nodes?
                //first find the op node corresponding to this mrrg node
                OpGraphOp *second_op = getOpNodePtr(opgraph, fu);

                //keep track of op's MRRG Node
                MRRGNode *first_MRRGNode = getMappedMRRGNode(op);

                //swap these two nodes
                OpMapping oldmap_x = ripUpOp(op);
                OpMapping oldmap_y = ripUpOp(second_op);

                bool success_x = placeOp(op, fu);
                bool success_y = placeOp(second_op, first_MRRGNode);

                //route
                if (success_x)
                {
                    bool routed = routeOp(op, mrrg);
                    assert(routed);
                }

                if (success_y)
                {
                    bool routed = routeOp(second_op, mrrg);
                    assert(routed);
                }

                float new_cost = getCost(mrrg);
                float delta_cost = new_cost - old_cost;
                if (accept(delta_cost, temperature))
                {
                    total_accepted++;
                }
                else //restore the opgraph
                {
                    ripUpOp(op);
                    ripUpOp(second_op);
                    restoreOp(oldmap_x);
                    restoreOp(oldmap_y);
                }
            }
        }
    }

    

    bool mrrg_overuse = checkOveruse(mrrg);
    float opgraph_cost = getCost(opgraph);

     LOG(DUMPMAPPING)<<printMapping()<<"\n";

    if (!mrrg_overuse)
    {
        cout << "MRRG OVERUSED!" << endl;
    }
    else
    {
        cout << "MRRG NOT OVERUSED!" << endl;
    }

    if (mrrg_overuse && opgraph_cost < INFINITY)
    {
        cout << "mrrg cost: " << getCost(mrrg) << endl;
        cout << "mrrg size: " << mrrg->routing_nodes.size() << endl;
        cout << "temp: " << temperature << endl;
        cout << "pfactor: " << pfactor << endl;
        *accept_rate = (float)total_accepted / (float)total_tries;
        return true;
    }
    *accept_rate = (float)total_accepted / (float)total_tries;
    return false;
}


 std::vector<std::pair<MRRGNode *, int>> LISAMapper::candidate_fu_intersect(std::vector<std::pair<MRRGNode *, int>> v1, std::vector<std::pair<MRRGNode *, int>> v2)
    {
        std::vector<std::pair<MRRGNode *, int>> result;

        for (auto &p : v1)
        {
            for (auto &n : v2)
            {
                if (p.first == n.first)
                    result.push_back({p.first, p.second + n.second});
            }
        }

        return result;
    }

    // This topological sort algorithm uses a depth first search
    // NB: The graph MUST be acyclic, else explosions may happen...
    void LISAMapper::topological_sort(OpGraph *g)
    {
        std::vector<OpGraphOp *> L;

        // Unmark all nodes
        std::map<OpGraphOp *, bool> mark;
        for (auto &n : g->op_nodes)
        {
            mark[n] = false;
        }

        OpGraphOp *n;
        while ((n = topological_sort_check_marked(g->op_nodes, &mark)))
        {
            topological_sort_visit(&L, &mark, n);
        }

        g->op_nodes = L;
    }
    /*
bool test_compare(float penalty_factor, MRRGNode* a, MRRGNode* b)
{
    return getCost(a) < getCost(b);
}

bool compare_op_cost(float penalty_factor, OpGraphOp* a, OpGraphOp* b)
{
    return a->getCost(penalty_factor) > b->getCost(penalty_factor);
}
*/
    float LISAMapper::updateTemperature(float t, float acceptance_rate)
    {
        if (acceptance_rate > 0.96)
        {
            return t * 0.5;
        }
        else if (acceptance_rate > 0.8)
        {
            return t * 0.9;
        }
        else if (acceptance_rate > 0.15)
        {
            return t * 0.98;
        }
        else
        {
            return t * 0.95;
        }
    }

    // generate a random FU
    MRRGNode *  LISAMapper::getRandomFU(MRRG *mrrg, OpGraphOp *op)
    {
        std::vector<MRRGNode *> candidates;
        for (auto &fu : mrrg->function_nodes)
        {
            if (fu->canMapOp(op) && occupancy[fu] == 0)
            {
                candidates.push_back(fu);
            }
        }

        if (candidates.size() < 1)
        {
            cout << "Could not place: " << *op << endl;
            assert(candidates.size() > 0);
        }

        std::random_device r;
        std::default_random_engine e1(r());
        std::uniform_int_distribution<int> uniform_dist(0, candidates.size()-1);
        int mean = uniform_dist(e1);

        return candidates[mean];
    }

    MRRGNode *  LISAMapper::getCloseRandomFU(MRRG *mrrg, OpGraphOp *op, MRRGNode * old_fu, int max_physical_dis ,  int max_temp_dis )
    {
        std::vector<MRRGNode *> opt_candidates;
        std::vector<MRRGNode *> all_candidates;

        int old_cycle = old_fu->cycle;
        std::set<int> allowed_cycle ;
        for(int start_cycle = old_cycle - max_temp_dis ; start_cycle <= old_cycle + max_temp_dis;start_cycle ++ ){
            allowed_cycle.insert( (start_cycle + II_)%II_);
        }

        LOG(RANDOMFU)<<"select CloseRandomFU for op"<<op->name<<" old fu"<<old_fu->getFullName();
        // assert(false);
        for (auto &fu : mrrg->function_nodes)
        {
            mpos mmmm{fu->x_, fu->y_, fu->cycle};
            bool exist = false;
            for(auto m :  mapped_pos){
                if(m.x == mmmm.x && mmmm.y == m.y  && mmmm.t == m.t){
                    exist = true;
                    break;
                }
            }
            if (fu->canMapOp(op) && occupancy[fu] == 0 )
            {
                all_candidates.push_back(fu);
                if(  std::abs(old_fu->x_ - fu->x_) + std::abs(old_fu->y_ - fu->y_) <= max_physical_dis  
                 && allowed_cycle.find(fu->cycle)!= allowed_cycle.end()){
                     opt_candidates.push_back(fu);
                 }
            }
        }
        if (all_candidates.size() < 1)
        {
            cout << "Could not place: " << *op << endl;
            assert(all_candidates.size() > 0);
        }
        std::random_device r;
        std::default_random_engine e1(r());
       
        if(opt_candidates.size() == 0){
            auto cal_dist = [](MRRGNode * old_fu, MRRGNode * new_fu){
                return std::abs(old_fu->x_ - new_fu->x_) + std::abs(old_fu->y_ - new_fu->y_) + std::abs(int(old_fu->cycle) - int(new_fu->cycle) );
            };
            std::sort(all_candidates.begin(), all_candidates.end(), [&](MRRGNode * a, MRRGNode * b){
                return cal_dist(old_fu, a) < cal_dist(old_fu, b);
            });
            std::uniform_int_distribution<int> uniform_dist(0, all_candidates.size()-1);
            int mean = uniform_dist(e1);
            auto sel_fu = all_candidates[mean];
            LOG(RANDOMFU)<<"  select fu:"<<sel_fu->getFullName();
            return sel_fu;
        }
        
        std::uniform_int_distribution<int> uniform_dist(0, opt_candidates.size()-1);
        int mean = uniform_dist(e1);
        auto sel_fu = opt_candidates[mean];
        LOG(RANDOMFU)<<"  select fu:"<<sel_fu->getFullName();
        return sel_fu;
        
    }
   

    template <typename type>
    std::vector<type> LISAMapper::vector_intersect(std::vector<type> v1, std::vector<type> v2)
    {
        std::sort(v1.begin(), v1.end());
        std::sort(v2.begin(), v2.end());

        std::vector<type> v_intersection;

        std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(v_intersection));

        return v_intersection;
    }

    std::vector<MRRGNode *> LISAMapper::filterCanMapOp(std::vector<MRRGNode *> fus, OpGraphOp const *op)
    {
        std::vector<MRRGNode *> result;
        for (auto &f : fus)
        {
            if (f->canMapOp(op))
            {
                result.push_back(f);
            }
        }

        return result;
    }



    // from herr, that is the new thingd that lisa added.

    std::map<int, pos3d> LISAMapper::dumpMapping(std::shared_ptr<OpGraph> opgraph,  int &max_latency){

        std::map<int, pos3d> dumpmapping;
        max_latency = 0;
        
       
        
        auto sorted_nodes = opgraph->op_nodes;
        std::sort (sorted_nodes.begin(), sorted_nodes.end(), 
        [&, this](OpGraphOp* a, OpGraphOp* b) {return node_asap_ [node_to_id_[a]] < node_asap_ [node_to_id_[b]]; 
        });
        LOG(CMAP)<<"*******dump mapping\n";
        for(auto node: sorted_nodes){
            if(mapping.find(node) == mapping.end() || mapping[node].size() == 0){
                continue;
            }
            auto m = mapping[node];
             LOG(CMAP)<<"node name:"<<node->name;
            if(dynamic_cast<OpGraphVal*>(node)) continue;
            auto node_op = dynamic_cast<OpGraphOp*>(node);
            assert(node_op);
            auto mrrg_node = m.front();
            LOG(CMAP)<<" mrrg_node:"<<mrrg_node->x_<<" "<<mrrg_node->y_<<" "<<mrrg_node->cycle;
            int input_max_lat = 0;
            for(auto val:  node_op->input){
                if(!val) continue; // remove backtrack edge, as the cgra-me does not support keeping data into register.
                int index = 0;
                bool find_input = false;
                for(auto output: val->output){
                   if(output == node ){
                    //    std::vector<OpGraphOp*> trans_input_op ;
                    //    for(auto trans_input_value : output->input){
                    //        if(!trans_input_value) continue;
                    //        trans_input_op.push_back(trans_input_value->input);
                    //    }
                    //    if(std::find(trans_input_op.begin(), trans_input_op.end(), node_op)==trans_input_op.end()) {
                            find_input = true;
                            break;
                    //    }
                       
                   }
                   index ++;
                }

                // for(int i = 0; i <  val->output.size(); i++){
                //     LOG(CMAP)<<"output node "<<i<<" "<<(val->output)[i]->name<< " latency: "<<(val->output_latency)[i];
                // }
                LOG(CMAP)<<"( input val "<<val->name<<" output size"<<val->output.size()<<" output latency size"<<val->output_latency.size();
                auto input_op = val->input;
                LOG(CMAP)<<" input_op name:"<<input_op->name<<" ";
                LOG(CMAP)<<dumpmapping[node_to_id_[input_op]].toStr();
                int lat = dumpmapping[node_to_id_[input_op]].t;
                LOG(CMAP)<<" lat:"<<lat<<")\n";
                input_max_lat = std::max(lat, input_max_lat);
                if(find_input  ){
                    
                    if(val->output_latency.size() == 0) continue;
                    if(dumpmapping.find(node_to_id_[input_op]) == dumpmapping.end()) continue;
                    LOG(CMAP)<<" output_latency:"<<(val->output_latency)[index]<<"\n";

                    int route_lat = (val->output_latency)[index];
                     //cgra-me did not add latency for function unit. Below code is for this bug.
                     auto op_code = val->input->opcode;
                    //  if(lat_opcode.find(op_code) != lat_opcode.end()){
                    //      route_lat += 1;
                    //  }
                     lat += route_lat;
                }
                LOG(CMAP)<<" lat:"<<lat<<")\n";
                input_max_lat = std::max(lat, input_max_lat);
            }
            // FIXME: is this right?
            int cycle = int((input_max_lat+1)/II_);
            int final_lat = cycle * II_ + mrrg_node->cycle;
            while(final_lat < input_max_lat){
                final_lat += II_; // due to cgra_me design bugs (cannot get output latency), it has a very low chance to get wrong latecny.
            }
            LOG(CMAP)<<" max lat"<<input_max_lat<<" final_lat:"<<final_lat<<"\n";
            dumpmapping[node_to_id_[node_op]]=pos3d{mrrg_node->x_, mrrg_node->y_, final_lat};
            LOG(CMAP)<<"final pos"<<dumpmapping[node_to_id_[node_op]].toStr()<<"\n";
        }
        // std::cout<<"finish";

        // max latency
        for(auto m: dumpmapping){
            int lat = m.second.t;
            max_latency = lat >  max_latency? lat: max_latency;
                // std::cout<<" value:"<<val->name<<" latency"<<lat<<std::endl;
        }
        max_latency = max_latency + 1;
        // std::cout<<"max latency"<<max_latency<<std::endl;


     
        // std::cout<<"dumpmapping2 "<<dumpmapping.size()<<"\n";
        return dumpmapping;
    }   

    MRRGNode *  LISAMapper::getLISAFU(std::shared_ptr<OpGraph> opgraph,MRRG *mrrg, OpGraphOp *op, int accpepted, int total_tried, int num_swap){
        auto & op_label = dfg_label_->at(node_to_id_[op]);

        std::set<OpGraphOp*> mapped_node;
        int temp_max_lat;// not usefule, just for call dumpMapping;
        auto dumped_mapping = dumpMapping(opgraph, temp_max_lat);
        // we should calculate based on mapped node
        for(auto node: node_to_id_){
            if(mapping[node.first].size()!= 0){
                mapped_node.insert(node.first);
            }
        }

        if(mapped_node.size() == 0){
            //should select fu that II is 1
            auto candidates = getDesiredFu(mrrg, 0 , 0, op);
            std::random_device r;
            std::default_random_engine e1(r());
            std::uniform_int_distribution<int> uniform_dist(0, candidates.size()-1);
            int mean = uniform_dist(e1);

            return candidates[mean];
        }       

        std::vector<MRRGNode *> candidates;
        for (auto &fu : mrrg->function_nodes)
        {
            mpos mmmm{fu->x_, fu->y_, fu->cycle};
            bool exist = false;
            for(auto m :  mapped_pos){
                if(m.x == mmmm.x && mmmm.y == m.y  && mmmm.t == m.t){
                    exist = true;
                    break;
                }
            }
            if (fu->canMapOp(op) && occupancy[fu] == 0 )
            {
                candidates.push_back(fu);
            }
        }

        // How to select FUs according to labels is a critical question. 
        //Give up this one: Let us use schedule order to filter some FUs, then we can use communication and association to select FU.
        // Trying make sure the node is scheduled as soon as possible, as long it satisifes schedule order constraints.
        // TODO: make every node as soon as  possible might be not accurate, check how to improve this. Maybe we can check data routing somehow.
        // 
        
        int schedule_order =   op_label.schedule_order;
        auto interval = getIntervalByScheduleOrder(opgraph, dumped_mapping, op, schedule_order);
         // if we only use label in initial mapping, then the second value is not useful
        int early_II =  interval.first;
        int start_II = 0;
        if(early_II +1 == II_) {
            start_II = 0;
        }else{
            start_II = early_II +1;
        }
        //  this sort is not used
        // std::sort(candidates.begin(), candidates.end(),[&](MRRGNode * a, MRRGNode * b){
        //     int a_t= a->cycle, b_t = b->cycle;
        //     if(a_t >= start_II && b_t >= start_II){
        //         return a_t < b_t;
        //     }else if(a_t <= start_II && b_t <= start_II){
        //         return a_t < b_t;
        //     }else{
        //         return a_t > b_t;
        //     }
        // });
        // std::cout<<"candidats size "<<candidates.size()<<"\n";
       
        assert(candidates.size() > 0);
        std::map<MRRGNode *, int> comm_cost;
        // if(!finish_init){
        //     getCostByComm(mrrg,opgraph, dumped_mapping, candidates, op );
        // }else{
        for(auto candi: candidates){
            comm_cost.emplace(candi,0);
        }
        // }
        
        auto samelevel_node_cost=  getCostForSameLevelNode(opgraph, dumped_mapping, candidates, op );
        auto ass_cost = getCostByAssociation(opgraph, dumped_mapping, candidates, op, start_II);
        
        std::stringstream output;

         //timing cost: this is to evaluate that whether we can route data in next iteration(for cgra_me only, 
         //for other frameworks, we need to add  simlar function but specfic to framework).
         // if not, we need to add cost for this.
        std::map<MRRGNode*, int> timing_cost ;
        auto & ass = dfg_label_->at(node_to_id_[op]).association;
        for(auto node: candidates){
            int ii_value = node->cycle;
            int t_cost = 0;
            for(auto node_ass: ass){
                int node_id = node_ass.first;
                if(dumped_mapping.find(node_id) == dumped_mapping.end()) continue;
                auto & m = dumped_mapping[node_id];
                int spatial_diff =  (std::abs(node->x_ - m.x) + std::abs(node->y_ - m.y));
                int tempoal_diff = MOD(m.t - ii_value);
                if(spatial_diff > tempoal_diff){
                    t_cost += spatial_diff - tempoal_diff ;
                }
                
            }
            timing_cost.emplace(node, t_cost);
        }
        output<<"early II: "<<early_II<<"\n";
        output<<"start II: "<<start_II<<"\n";
        output<<"\n\ttiming cost:";
        for(auto node: timing_cost){
            output<<"("<<node.first->getFullName()<<","<<node.second<<") ";
        }

        output<<"\n\tassociation cost:";
        for(auto node: ass_cost){
            output<<"("<<node.first->getFullName()<<","<<node.second<<") ";
        }
        output<<"\n\tcommunication cost:";
        for(auto node: comm_cost){
            output<<"("<<node.first->getFullName()<<","<<node.second<<") ";
        }

         output<<"\n\tstart_node_cost cost:";
        for(auto node: samelevel_node_cost){
            output<<"("<<node.first->getFullName()<<","<<node.second<<") ";
        }
        
       


        std::map<MRRGNode *, int> node_cost ;
        for(auto node : candidates){
            int cost = timing_cost[node]  + ass_cost[node] + samelevel_node_cost[node] ;
            node_cost.emplace(node, cost);
        }

        output<<"\n\t total cost:";
        for(auto node: node_cost){
            output<<"("<<node.first->getFullName()<<","<<node.second<<") ";
        }
        
        std::sort(candidates.begin(), candidates.end(), [&](MRRGNode * a, MRRGNode * b){
            return node_cost[a] < node_cost[b];
        });

        
        int min_cost = node_cost[candidates.front()];
        int max_cost = node_cost[candidates.back()];
        int max_diff = max_cost -  min_cost;

        std::random_device rd{};
        std::mt19937 gen{rd()};
        double deviation = 1;
        // std::cout<<"1,2,3:"<<accpepted<<","<<total_tried<<","<<num_swap<<"\n";
        if(total_tried > 0){
            deviation = 0.1 * total_tried - accpepted;
            deviation = std::max(deviation, 1.0);
                
        }
       
        std::normal_distribution<> d{0,deviation};
        int selected_cost = std::abs(std::round(d(gen)));
        if(selected_cost > max_diff)  selected_cost = max_diff;
        selected_cost  = min_cost + selected_cost;

        if(!finish_init) selected_cost = min_cost;

        // 
        std::vector<MRRGNode *> suitable_candidates;
        while(suitable_candidates.size() == 0){
            for(auto c : node_cost){
                if(c.second == selected_cost){
                    suitable_candidates.push_back(c.first);
                }
            }
            selected_cost --;
        }
        
        // std::cout<<"min_candidates size:"<<min_candidates.size()<<"\n";
        std::random_device r;
        std::default_random_engine e1(r());
        std::uniform_int_distribution<int> uniform_dist(0, suitable_candidates.size()-1);
        int mean = uniform_dist(e1);
        auto selected_fu = suitable_candidates[mean];


        if(!is_training){
            LOG(LISAFU)<<"curr mapping: "<<printMapping();
            LOG(LISAFU)<<"label of "<<*op<<"  " <<op_label.toStr();
            LOG(DLABEL)<<output.str();
            LOG(LISAFU)<<"best mapping: "<<lisa_ctrl->bestMappingToStr();
            LOG(LISAFU)<<" min cost:"<< min_cost<<"   selected cost: "<<selected_cost +1 <<"\n";
            LOG(LISAFU)<<"lisa op:"<<op->name<<" selectFU: "<<selected_fu->getFullName()<<"\n";
        }

       
        
        return selected_fu;

    } 
    std::map<MRRGNode *, int> LISAMapper::getCostByComm(MRRG *mrrg, std::shared_ptr<OpGraph> opgraph, std::map<int, pos3d> & dumped_mapping, std::vector<MRRGNode *> candidates, OpGraphOp *op ){
        /*
        0,0     0,1     1,0     1,1
        #       #       #       #
        *       *       used    #
        #       candi   #       #
        *       used    *       #
        # and * means emprt one, while * means routing resource.
        the number of routing resource should be 4
        */
         
         auto find_routing_resource_ = [&](MRRGNode * a , bool up_directoin){
            int curr_cycle = a->cycle;
            std::set<MRRGNode * > overall_routing_resource;
            std::set<MRRGNode * > curr_routing_resource;
            curr_routing_resource.insert(a);
            while(true){
                std::set<MRRGNode * > new_node;
                if(curr_routing_resource.size() ==  0){
                    break;
                }

                //check whether blocked
                bool path_blocked = false;
                for(auto node: curr_routing_resource){
                    int x = node->x_, y = node->y_;
                    int blocked_side = 0;
                    if(x-1 < 0)  blocked_side ++;
                    if(x+1 >= cgra_x_)  blocked_side ++;
                    if(y-1 < 0)  blocked_side ++;
                    if(y+1 >= cgra_y_)  blocked_side ++;
                    if( blocked_side >= 2) {
                        break;
                        path_blocked = true;
                    }    

                }
                if(path_blocked){
                    break;
                }

                curr_cycle = MOD(curr_cycle);
                int next_cycle ;
                if(up_directoin) {next_cycle = curr_cycle - 1;}  else {next_cycle = curr_cycle + 1 ;}
                next_cycle = MOD(next_cycle);

                for(auto node: curr_routing_resource){
                    if( node->cycle == curr_cycle){
                        //find the nieghbors
                        int x = node->x_, y = node->y_;
                        auto temp_node = getRoutingNode(mrrg, x-1, y, next_cycle); if(temp_node && occupancy[temp_node] == 0) new_node.insert(temp_node);
                        temp_node = getRoutingNode(mrrg, x+1, y, next_cycle); if(temp_node && occupancy[temp_node] == 0) new_node.insert(temp_node);
                        temp_node = getRoutingNode(mrrg, x, y-1, next_cycle); if(temp_node && occupancy[temp_node] == 0) new_node.insert(temp_node);
                        temp_node = getRoutingNode(mrrg, x, y+1, next_cycle); if(temp_node && occupancy[temp_node] == 0) new_node.insert(temp_node);
                    }
                }
                curr_routing_resource.clear();
                for(auto node: new_node){
                    if(overall_routing_resource.find(node) == overall_routing_resource.end()){
                        curr_routing_resource.insert(node);
                        overall_routing_resource.insert(node);
                    }
                   
                }
                if(up_directoin) {curr_cycle --;}  else {curr_cycle ++ ;}
                
            }
            return overall_routing_resource;
        };

        std::map<MRRGNode *, int> candi_routing_resource_num;
        for(auto candi : candidates){
            auto up_routing_resource = find_routing_resource_(candi, true);
            auto down_routing_resource = find_routing_resource_(candi, false);
            std::set<MRRGNode *> overall_routing_resource;
            for(auto r :  up_routing_resource){
                overall_routing_resource.insert(r);
            }
            for(auto r :  down_routing_resource){
                overall_routing_resource.insert(r);
            }
            candi_routing_resource_num[candi] =  overall_routing_resource.size();
        }
        return candi_routing_resource_num;
    }



    MRRGNode* LISAMapper::getRoutingNode(MRRG *mrrg, int x, int y, int t){
        for(auto node: mrrg->routing_nodes){
            if(node->x_ == x && node->y_ == y  && node->cycle == t && node->latency == 1){
                return node;
            }
        }
        return NULL;
    }

    std::map<MRRGNode *, int>   LISAMapper::getCostByAssociation
    (std::shared_ptr<OpGraph> opgraph, std::map<int, pos3d> & dumped_mapping, std::vector<MRRGNode *> candidates, OpGraphOp *op ,  int start_II){
        auto & ass = dfg_label_->at(node_to_id_[op]).association;
        int earliest_execution_time = 0;
        int op_id = node_to_id_[op];
        for(auto e: lisa_edges_){
            if(e.second == op_id){
                int temp_time = dumped_mapping[e.first].t;
                earliest_execution_time = std::max(earliest_execution_time, temp_time);
            }
        }
        int earliest_execution_time_II = MOD(earliest_execution_time);

    
        auto get_cost = [&, this ](MRRGNode * a) {
            int mrrgnode_II =  a->cycle;
            int total_spatial_cost = 0, total_temp_cost = 0;
            int guess_time = 0;

            if(mrrgnode_II < earliest_execution_time_II){ 
                guess_time = earliest_execution_time + ( mrrgnode_II + II_) - earliest_execution_time_II;
            }else{
                guess_time = earliest_execution_time + mrrgnode_II + - earliest_execution_time_II;
            }
            int mapped_ass_node = 0;
            for(auto node_ass: ass){
                int node_id = node_ass.first;
                if(dumped_mapping.find(node_id) == dumped_mapping.end()) continue;
                mapped_ass_node++;
                auto & m = dumped_mapping[node_id];
                total_spatial_cost += std::abs(std::abs(node_ass.second.first) - (std::abs(a->x_ - m.x) + std::abs(a->y_ - m.y)));
                total_temp_cost += std::abs(m.t - guess_time);
                
            }
            if(mapped_ass_node == 0) return std::make_pair(0,0);
            return std::make_pair(int(total_spatial_cost/mapped_ass_node), int(total_temp_cost/mapped_ass_node));
        };

        std::map<MRRGNode *, int> candi_ass_cost;
        for(auto candi: candidates){
            auto cost = get_cost(candi);
            candi_ass_cost.emplace(candi, cost.first + cost.second);
        }

        return candi_ass_cost;
    }


    std::map<MRRGNode *, int> LISAMapper::getCostForSameLevelNode
    (std::shared_ptr<OpGraph> opgraph, std::map<int, pos3d> & dumped_mapping, std::vector<MRRGNode *> candidates, OpGraphOp *op ){
        int node_id = node_to_id_[op];
        std::map<MRRGNode *, int> node_cost;
         for(auto node: candidates){
            node_cost[node] = 0;
        }

        // if(! lisa_ctrl->isStartNode(node_id)) return node_cost;

        auto relevant_samv_level_nodes = lisa_ctrl->getSameLevelNodes(node_id);
        std::set<int> mapped_sameLevel_nodes;
        for(auto rele_node: relevant_samv_level_nodes){
            if(dumped_mapping.find(rele_node)!= dumped_mapping.end()){
                mapped_sameLevel_nodes.insert(rele_node);
            }
        }
        // if(mapped_sameLevel_nodes.size()== 0) assert(false);
        auto & dist_label = dfg_label_->at(node_id).sameLevel_node_distance; 
        auto cal_cost = [&](MRRGNode * a){
            if(mapped_sameLevel_nodes.size()==0) return 0;
            int total_cost = 0;
            for(auto node: mapped_sameLevel_nodes){
                auto & m = dumped_mapping[node];
                int dist = std::abs(a->x_ - m.x) + std::abs(a->y_ - m.y);
                total_cost += std::abs(dist_label[node] - dist);
            }
            return (int)(total_cost/(mapped_sameLevel_nodes.size()));
        };

        for(auto node: candidates){
            node_cost[node] = cal_cost(node);
        }

        return node_cost;
    }

    std::pair<int,int> LISAMapper::getIntervalByScheduleOrder
    (std::shared_ptr<OpGraph> opgraph, std::map<int, pos3d> & dumped_mapping, OpGraphOp *op,int scheduler_order ){
        // return value is on II.
        int start_time = -1, end_time = -1;

        auto  early_ops_comp =   [ &](OpGraphOp * a, OpGraphOp * b) {
            return dumped_mapping[node_to_id_ [a]].t  > dumped_mapping [node_to_id_ [b]].t;
        };

        auto  late_ops_comp =    [ &](OpGraphOp * a, OpGraphOp * b) {
            return dumped_mapping[node_to_id_ [a]].t < dumped_mapping [node_to_id_ [b]].t;
        };

        std::set<OpGraphOp *, decltype(early_ops_comp)> early_ops(early_ops_comp);
        std::set<OpGraphOp *, decltype(late_ops_comp)> late_ops(late_ops_comp);
        // std::cout<<"this node"<<scheduler_order<<"\n";
        for(auto node: opgraph->op_nodes){
            // std::cout<< node<<" order: "<<dfg_label_->at(node_to_id_[node]).schedule_order<<"\n"; 
            if(dfg_label_->at(node_to_id_[node]).schedule_order < scheduler_order  && mapping[node].size()!=0 ) early_ops.insert(node);
            if(dfg_label_->at(node_to_id_[node]).schedule_order > scheduler_order  && mapping[node].size()!=0) late_ops.insert(node);
        }
       
        if(early_ops.size() != 0){
           start_time = mapping[*(early_ops.begin())].front()->cycle;
        }
        if(late_ops.size() != 0){
           end_time = mapping[*(late_ops.begin())].front()->cycle;;
        }

        return std::make_pair(start_time, end_time);

    }

    std::vector<MRRGNode*> LISAMapper::getDesiredFu(MRRG *mrrg, int start_II, int end_II, OpGraphOp *op){
        std::vector<MRRGNode *> candidates;
        for (auto &fu : mrrg->function_nodes)
        {
            mpos mmmm{fu->x_, fu->y_, fu->cycle};
            bool exist = false;
            for(auto m :  mapped_pos){
                if(m.x == mmmm.x && mmmm.y == m.y  && mmmm.t == m.t){
                    exist = true;
                    break;
                }
            }
            if (fu->canMapOp(op) && occupancy[fu] == 0 && !exist)
            {
                if(fu->cycle >= start_II && fu->cycle <= end_II){
                    candidates.push_back(fu);
                }
            }
        }

        std::sort(candidates.begin(), candidates.end(),
        [](MRRGNode * a, MRRGNode * b) {
            return a->cycle< b->cycle || ( a->cycle == b->cycle && a->x_+ a->y_ < b->x_+ b->y_) 
            || ( a->cycle == b->cycle && a->x_+ a->y_ == b->x_+ b->y_ && a->x_ < b->x_); 
        });
        assert(candidates.size() > 0);
        return candidates;
    }
   