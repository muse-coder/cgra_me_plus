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
#include <CGRA/AnnealMapper.h>

using std::cout;
using std::endl;

#define ALLOW_MULTIPLE_PLACEMENT false/
#define SIMPLE_ANNEAL_SCHEDULE
#define CALCULATE_INITIAL_TEMPERATURE
// Debugging output flags
//#define ANNEAL_DEBUG
//#define PLACEMENT_DEBUG
// #define DEBUG_ROUTING

AnnealMapper::AnnealMapper(std::shared_ptr<CGRA> cgra, int timelimit, const std::map<std::string, std::string> &args)
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
AnnealMapper::AnnealMapper(std::shared_ptr<CGRA> cgra, int timelimit, int rand_seed, float initial_pfactor, float pfactor_factor, float const_temp_factor, int swap_factor, float cold_accept_rate)
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

AnnealMapper::AnnealMapper(std::shared_ptr<CGRA> cgra)
    : AnnealMapper(cgra, 0, 0, 0.001, 1.05, 0.99, 100, 0.01)
{
}

float AnnealMapper::getCost(MRRGNode *n)
{
    float base_cost;

    if (n->type == MRRG_NODE_FUNCTION)
        base_cost = 2.0;
    else
        base_cost = 1.0;

    return base_cost * occupancy[n] + (occupancy[n] <= n->capacity ? 0.0 : (occupancy[n] - n->capacity) * pfactor);
}

float AnnealMapper::getCost(MRRG *mrrg)
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

float AnnealMapper::getCost(OpGraph *opgraph)
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

float AnnealMapper::getCost(OpGraphNode *n)
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

bool AnnealMapper::checkOveruse(MRRG *mrrg)
{
    bool result = true;
    for (auto &node : mrrg->function_nodes)
    {
        if (occupancy[node] > node->capacity)
        {
            LOG(OVERUSE) << *node << " is overused. (" << occupancy[node] << "/" << node->capacity << ")\n";
            result = false;
        }
    }

    for (auto &node : mrrg->routing_nodes)
    {
        if (occupancy[node] > node->capacity)
        {
            LOG(OVERUSE) << *node << " is overused. (" << occupancy[node] << "/" << node->capacity << ")\n";
            result = false;
        }
    }

    return result;
}

OpGraphOp *AnnealMapper::getOpNodePtr(OpGraph *opgraph, MRRGNode *n)
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

void AnnealMapper::mapMRRGNode(OpGraphNode *opnode, MRRGNode *n)
{
    try
    {
        mapping[opnode].push_back(n);
    }
    catch (const std::exception &e)
    {
        throw cgrame_error(std::string("AnnealMapper Exception Thrown by: [") + e.what() + "] at File: " + std::string(__FILE__) + " Line: " + std::to_string(__LINE__));
    }
    occupancy[n]++;
}

void AnnealMapper::mapAllMRRGNodes(OpGraphNode *opnode, std::vector<MRRGNode *> nodes)
{
    // map all nodes
    for (auto &n : nodes)
    {
        mapMRRGNode(opnode, n);
    }
}

void AnnealMapper::unmapMRRGNode(OpGraphNode *opnode, MRRGNode *n)
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
            occupancy[n]--;
            assert(occupancy[n] >= 0);
            mapping[opnode].erase(iter);
        }
    }
    catch (const std::exception &e)
    {
        throw cgrame_error(std::string("AnnealMapper Exception Thrown by: [") + e.what() + "] at File: " + std::string(__FILE__) + " Line: " + std::to_string(__LINE__));
    }
}

std::vector<MRRGNode *> AnnealMapper::unmapAllMRRGNodes(OpGraphNode *opnode)
{
    std::vector<MRRGNode *> result;

    // save mapping
    result = mapping[opnode];

    // unmap all nodes
    for (auto &n : mapping[opnode])
    {
        occupancy[n]--;
        assert(occupancy[n] >= 0);
    }

    mapping[opnode].clear();

    return result;
}

/**
  Returns the MRRGNode that the Op is mapped to, NULL if unmapped
 **/
MRRGNode *AnnealMapper::getMappedMRRGNode(OpGraphOp *op)
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
bool AnnealMapper::routeVal(OpGraphVal *val)
{

    assert(val);
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
    // assert(false);
    return true;

}

float AnnealMapper::getTotalOpCost(OpGraphOp *op)
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
OpMapping AnnealMapper::ripUpOp(OpGraphOp *op, float *cost)
{
    OpMapping result;

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

void AnnealMapper::restoreOp(OpMapping oldmap)
{
    for (auto m = oldmap.mapping.begin(); m != oldmap.mapping.end(); ++m)
    {
        mapAllMRRGNodes(m->first, m->second);
    }
}

// get a random unoccupied  FU
MRRGNode *AnnealMapper::getRandomUnoccupiedFU(MRRG *mrrg, OpGraphOp *op)
{
    int min_cycle = 0;
    for(auto input_val: op->input){
        auto op =  input_val->input;
         if (mapping[op].size() != 0){
            int cyclye = mapping[op].front()->cycle;
            min_cycle =  cyclye>min_cycle? cyclye:min_cycle;
         }

    }
    std::vector<MRRGNode *> candidates;
    std::vector<MRRGNode *> badcandidates;
    for (auto &fu : mrrg->function_nodes)
    {
        if (fu->canMapOp(op) && occupancy[fu] == 0)
        {
            if(fu->cycle < min_cycle){
                candidates.push_back(fu);
            }
            badcandidates.push_back(fu);
        }
    }
    if (candidates.size() < 1){
        candidates =  badcandidates;
    }else {
            return candidates[rand() % candidates.size()];

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
bool AnnealMapper::placeOp(OpGraphOp *op, MRRGNode *n)
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

inline float AnnealMapper::updateTempConst(float t)
{
    return t * const_temp_factor; //0.999;
}

bool AnnealMapper::routeOp(OpGraphOp *op, MRRG *mrrg)
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
    // LOG(SA)<<" route "<< op->name << " "<<result<<"\n";
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
Mapping AnnealMapper::mapOpGraph(std::shared_ptr<OpGraph> opgraph, int II, std::string arch_model_name)
{
    // get the mrrg object
    MRRG *mrrg = cgra->getMRRG(II).get();
//
//     for(int i = 0; i< systolic_array_x; i++){
//        for(int j=0; j< systolic_array_y;j++){
//            systolic_pe_index.emplace(std::make_pair(i,j),sys_arr[i][j]);
//        }
//    }
        arch_model_name_ = arch_model_name;


    this->this_mrr = mrrg;
//    if (arch_model_name_.find("systolic") != std::string::npos) {
//        this_mrr->makeSystolicArray(systolic_pe_index, systolic_array_x);
//    }
//
//     if (arch_model_name_.find("leftmostmemory") != std::string::npos) {
//        this_mrr->makeLeftMostMemoryAccess();
//    }

    MAX_OVERUSE = opgraph->op_nodes.size() * 100;

//    ***************initial MRRG and graph node***************
    // Create result obj
    Mapping mapping_result(cgra, II, opgraph);

    // Set the random seed
    srand(this->rand_seed);

     for(auto node: mrrg->routing_nodes){
        occupancy[node] = 0;
    }
    for(auto node: mrrg->function_nodes){
        occupancy[node] = 0;
    }

#ifdef ANNEAL_DEBUG
    ofstream anneal_debug;
    anneal_debug.open("anneal_debug.csv");
#endif

    // Sort the opgraph nodes
    // TODO: This does not work on graphs with back edges
    // topological_sort(opgraph);

    LOG(SA) << "Initial placement:" << "\n";

     std::vector<OpGraphOp*>  sorted_nodes;
    for(auto node: opgraph->op_nodes){
        if(find(sorted_nodes.begin(), sorted_nodes.end(),node) == sorted_nodes.end()){
            sorted_nodes.push_back(node);
        }
    }
    // randomized initial placement
    std::map<OpGraphNode*, MRRGNode*>  bestPlacement;
    int min_cost  = 10000000;
     bool routed = true;
//     搜索一百次找到当前最佳的placement
    for(int i  = 0 ; i < 100 ; i++){
        routed = true;
        std::map<OpGraphNode*, MRRGNode*>  currPlacement;

         for (auto &op : sorted_nodes)
        {
    #ifdef ALLOW_MULTIPLE_PLACEMENT
            MRRGNode *fu = getRandomFU(mrrg, op);
    #else
            MRRGNode *fu = getRandomUnoccupiedFU(mrrg, op);
    #endif
            LOG(SA) <<"iter:"<< i <<" " << *op << " : " << *fu << "\n";
            bool placed = placeOp(op, fu);
            currPlacement.emplace(op, fu);
            assert(placed);
        }

        for (auto &op :sorted_nodes)
        {
            routed = routeOp(op, mrrg);
            
            if(!routed){
//                std::cout<<"unrouted\n";
                break;
            }
        }

        if(routed) {
            int curr_cost = getCost(mrrg);
            if(curr_cost < min_cost){
                bestPlacement =  currPlacement;
                min_cost = curr_cost;
            }
        }

        for (auto &op :sorted_nodes)
        {
            ripUpOp(op);
        }

        if( i == 50) {
            allow_overuse =  true;//50次后再允许overuse
            LOG(SA) <<"iter:"<< i <<" set overuse\n" ;
        }
        
    }
    
    if(!routed) return mapping_result;
    for (auto &op : sorted_nodes)
    {
        ripUpOp(op);
    }
        for (auto &op : sorted_nodes)
    {
        placeOp(op, bestPlacement[op]);
    }
    for (auto &op :sorted_nodes)
    {
        routed = routeOp(op, mrrg);
        
        assert(routed);
    }

    allow_overuse =  false;

    // Verify initial place and route
    for (auto &op : sorted_nodes)
    {
        if (mapping[op].size() == 0)
        {
            cout << "No initial placement for: " << op->name << endl;
            assert(0);
        }
    }

    for (auto &val : sorted_nodes)
    {
        if (mapping[val].size() == 0)
        {
            cout << "No initial routing for: " << val->name << endl;
            assert(0);
        }
    }

#ifdef CALCULATE_INITIAL_TEMPERATURE
    /************************* TRY TO DO 100 iterations before annealing********/

    cout << "Finding delta Costs:" << endl
         << "INITIAL:" << endl
         << "OpGraph Cost: " << getCost(opgraph.get()) << endl;
    cout << "MRRG Cost: " << getCost(mrrg) << endl;
    cout << "Penalty Factor: " << pfactor << endl;

    float max_delta_cost = 0;
    float old_cost = getCost(mrrg);

    //do 100x of perturbation
    for (int i = 0; i < 100; i++)
    {
        //first get a random index
        int vector_size = opgraph->op_nodes.size();
        int index = rand() % vector_size;

        //perturb at this index
        OpGraphOp *op = opgraph->op_nodes[index];

        //substitute a new fu
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
                if(!success){
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
#ifndef ALLOW_MULTIPLE_PLACEMENT
                assert(fu->occupancy == 1);
#endif
                //swap, rip off two nodes?
                //first find the op node corresponding to this mrrg node
                OpGraphOp *second_op = getOpNodePtr(opgraph.get(), fu);

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

                //get new cost
                float new_cost = getCost(mrrg);
                float change_cost = abs(new_cost - old_cost);
                if (change_cost > max_delta_cost)
                    max_delta_cost = change_cost;

                //restore the opgraph
                ripUpOp(op);
                ripUpOp(second_op);
                restoreOp(oldmap_x);
                restoreOp(oldmap_y);
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
        float accept_rate = 0.0;
        float previous_cost = current_cost;
        LOG(SA) << "Annealing at:"
                << "\ttemp: " << temperature
                << "\tpfactor: " << pfactor;
        if (inner_place_and_route_loop(opgraph.get(), mrrg, temperature, &accept_rate))
        {
            mapping_result.setMapping(mapping);
            LOG(SA) << "MappingTime: " << (int)(getcurrenttime() - start_time);
            LOG(SA) << "MapperTimeout: 0";
            LOG(SA) << "Mapped: 1";
            mapping_result.setMapped(true);
            return mapping_result;
        }
        current_cost = getCost(mrrg);

#ifdef ANNEAL_DEBUG
        anneal_debug << temperature << ",";
        anneal_debug << current_cost << ",";
        anneal_debug << pfactor << ",";
        anneal_debug << endl;
#endif
        // TODO: Changed from 0.01
        //if(temperature  < 0.001 * mrrg->getCost(pfactor) / mrrg->routing_nodes.size())
        if (accept_rate < cold_accept_rate && current_cost >= previous_cost)
        {
#ifdef ANNEAL_DEBUG
            anneal_debug.close();
#endif
            LOG(SA) << "Mapper is Cold and no valid mapping was found.";
            LOG(SA) << "Current Temperature acceptance rate was: " << accept_rate;
            LOG(SA) << "Cold acceptance rate is: " << cold_accept_rate;
            LOG(SA) << "Current Cost is: " << current_cost;
            LOG(SA) << "Previous  Cost was: " << previous_cost;
            LOG(SA) << "Mapped: 0";
            return mapping_result;
        }
        LOG(SA) << "mrrg cost: " << getCost(mrrg);
        LOG(SA) << "mrrg size: " << mrrg->routing_nodes.size();
        LOG(SA) << "update temp. & pfactor";
        // update temperature
#ifdef SIMPLE_ANNEAL_SCHEDULE
        temperature = updateTempConst(temperature);
#else
        temperature = updateTemperature(temperature, (float)total_accepted / total_tries);
#endif
        // update overuse penalty
        pfactor = pfactor * pfactor_factor;

        current_time = getcurrenttime();

        LOG(SA) << "current run time: " << (int)(current_time - start_time);
    }

    cout << "MapperTimeout: 1" << endl;
    cout << "Mapped: 0" << endl;

    return mapping_result;
}

bool AnnealMapper::inner_place_and_route_loop(OpGraph *opgraph, MRRG *mrrg, float temperature, float *accept_rate)
{
    int num_swaps =  swap_factor  ;
    int total_accepted = 0;
    int total_tries = 0;

    for (int i = 0; i < num_swaps; i++)
    {
        // Get an op
        OpGraphOp *op = opgraph->op_nodes[rand() % opgraph->op_nodes.size()];

        // Get an fu
        MRRGNode *fu = getRandomFU(mrrg, op);
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
                    assert(success);
                }

                if (success)
                {
                    success = routeOp(op, mrrg);
                    if(!success){
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
