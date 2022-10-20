#include "netflow.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <map>

#define TOO_FEW_VERTICES "Too few vertices."
#define TOO_FEW_EDGES "Too few edges."
#define EDGE_WEIGHT_ZERO "Detected edge weight of 0."
#define EDGE_BAD_ENDPOINT "Edge interacts with nonexistent vertex."
#define SELF_LOOP "At least one self-loop."
#define MULTI_EDGES "Detected multi-edges."
#define NOT_ONE_SRC "Zero or more than one source."
#define NOT_ONE_SINK "Zero or more than one sink."

static std::vector<std::vector<int>> adjMatrix;
static std::vector<std::vector<int>> flow;
static std::vector<std::vector<int>> residual;
static unsigned source = 0 ;
static unsigned sink = 0;

static std::vector<Edge> input_vector;
static std::map<std::string, unsigned> match;
static std::map<unsigned, std::string> match_reversed;

void find_source_sink(const std::vector<Edge>& flowNetwork, unsigned numVertices)
{
    unsigned src_count = 0;
    unsigned sink_count = 0;
    source = numVertices ;
    sink = numVertices;
    
    for (unsigned i = 0; i < flowNetwork.size(); i++)
    {
        //std::cout << "from: " << flowNetwork.at(i).from << std::endl;
        unsigned src_itr = 0;
        unsigned sink_itr = 0;
        for (unsigned j = 0; j < flowNetwork.size(); j++)
        {
            //std::cout << "to: " << flowNetwork.at(j).to << std::endl;
            if (flowNetwork.at(i).from != flowNetwork.at(j).to)
            {
                src_itr++;
            }
            if (flowNetwork.at(i).to != flowNetwork.at(j).from)
            {
                sink_itr++;
            }
        }
        if (src_itr == flowNetwork.size())
        {
            if (flowNetwork.at(i).from != source)
            {
                source = flowNetwork.at(i).from;
                src_count++;
            }
        }
        if (sink_itr == flowNetwork.size())
        {
            if (flowNetwork.at(i).to != sink)
            {
                sink = flowNetwork.at(i).to;
                sink_count++;
            }
        }
    }
    if (src_count != 1)
    {
        throw std::runtime_error(NOT_ONE_SRC);
    }
    if (sink_count != 1)
    {
        throw std::runtime_error(NOT_ONE_SINK);
    }
}

void set_adjMatrix(const std::vector<Edge>& flowNetwork, unsigned numVertices)
{
    //initialize the adjacency matrix to be 0
    for(unsigned i = 0; i < numVertices; i++)
    {
        std::vector<int> from;
        adjMatrix.push_back(from);
        flow.push_back(from);
        for(unsigned j = 0; j < numVertices; j++)
        {
            adjMatrix.at(i).push_back(0);
            flow.at(i).push_back(0);
        }
    }

    for (unsigned i = 0; i < flowNetwork.size(); i++)
    {
        if (adjMatrix[flowNetwork.at(i).from][flowNetwork.at(i).to] != 0)
        {
            throw std::runtime_error(MULTI_EDGES);
        }
        adjMatrix[flowNetwork.at(i).from][flowNetwork.at(i).to] = flowNetwork.at(i).weight;
    }
    
}

bool BFS(unsigned numVertices, std::vector<int>& parent)
{
    //std::cout << "---- bfs starts -----" << std::endl;
    for (unsigned i = 0; i < numVertices; i++)
        {
            parent.at(i) = (-1);
        }
    //parent.at(source) = -1;
    std::vector<bool> discovered;
    for (unsigned i = 0; i < numVertices; i++)
    {
        discovered.push_back(false);
    }
    discovered.at(source) = true;

    std::queue<int> Queue;
    Queue.push(source);
    while (!Queue.empty())
    {
        unsigned curr = Queue.front();
        Queue.pop();
        for (unsigned i = 0; i < numVertices; i++)
        {
            //find neighbor
            if (residual.at(curr).at(i) != 0)
            {
                unsigned neighbor = i;
                if (discovered.at(neighbor) == false)
                {
                    discovered.at(neighbor) = true;
                    Queue.push(neighbor);
                    parent.at(neighbor) = curr;
                }
            }
        }
    }
    if (discovered.at(sink)==true)
    {
        return true;
    }
    return false;
}

unsigned find_bottle_neck(std::vector<int> parent)
{   
    unsigned temp = sink;
    unsigned previous = parent.at(sink);
    unsigned edge = residual.at(previous).at(sink);
    unsigned min = edge;
    while (temp != source)
    {
        previous = parent.at(temp);
        edge = residual.at(previous).at(temp);
        if (edge < min)
        {
            min = edge;
        }
        temp = previous;
    }
    return min;
}

void Ford_Fulkerson(unsigned numVertices)
{
    //create the residual graph
    for(unsigned i = 0; i < numVertices; i++)
    {
        std::vector<int> from;
        residual.push_back(from);
        for(unsigned j = 0; j < numVertices; j++)
        {
            residual.at(i).push_back(adjMatrix.at(i).at(j));
        }
    }
    //create an vector to store the path from BFS
    std::vector<int> parent;  //pointer?
    for (unsigned i = 0; i < numVertices; i++)
    {
        parent.push_back(-1);
    }
/*
    BFS(numVertices, parent);
    std::cout << "---- parent vector -----" << std::endl;
    
    for (unsigned i = 0; i < parent.size(); i++)
    {
        std::cout << i << " : " << parent.at(i) << std::endl;
    }
*/
    while (BFS(numVertices, parent))
    {
        // use parent[] to go through the path
        unsigned bottle_neck = find_bottle_neck(parent);
        //std::cout << "bottle_neck: " << bottle_neck << std::endl;
        unsigned temp = sink;
        while (temp != source)
        {
            unsigned previous = parent.at(temp);   
            //update the path flow
            //if it is a backward edge
            if (adjMatrix.at(previous).at(temp) == 0)
            {
                flow.at(temp).at(previous) -= bottle_neck;
                residual.at(previous).at(temp) -= bottle_neck;
                residual.at(previous).at(temp) += bottle_neck;
            }
            else
            {
                flow.at(previous).at(temp) += bottle_neck;
                residual.at(previous).at(temp) -= bottle_neck;
                residual.at(temp).at(previous) += bottle_neck;
            }
            temp = previous;
        }
    }
}

std::vector<Edge> solveNetworkFlow(
    const std::vector<Edge>& flowNetwork,
    unsigned numVertices)
{
    //input validation
    if (numVertices < 2)
    {
        throw std::runtime_error(TOO_FEW_VERTICES);
    }
    if (flowNetwork.empty())
    {
        throw std::runtime_error(TOO_FEW_EDGES);
    }

    for (unsigned i = 0; i < flowNetwork.size(); i++)
    {
        if (flowNetwork.at(i).weight == 0)
        {
            throw std::runtime_error(EDGE_WEIGHT_ZERO);
        }
        if (flowNetwork.at(i).from >= numVertices || flowNetwork.at(i).to >= numVertices)
        {
            throw std::runtime_error(EDGE_BAD_ENDPOINT);
        }
        if (flowNetwork.at(i).from == flowNetwork.at(i).to)
        {
            throw std::runtime_error(SELF_LOOP);
        }
    }
    set_adjMatrix(flowNetwork, numVertices);
    find_source_sink(flowNetwork, numVertices);
    Ford_Fulkerson(numVertices);
    std::vector<Edge> output;
    for (unsigned i = 0; i < numVertices; i++)
    {
        for (unsigned j = 0; j < numVertices; j++)
        {
            if(flow.at(i).at(j) != 0) 
            {
                output.push_back(Edge(i, j, flow.at(i).at(j)));
            }
        }
    }
    return output;
}



void assignCourses(
    std::vector<Instructor>& instructors,
    const std::vector<std::string>& courses)
{
    unsigned s = 0;
    unsigned num = instructors.size() + courses.size();
    unsigned t = num + 1;
    
    for (unsigned i = 0; i < instructors.size(); i++)
    {
        //std::cout << i << " : " << instructors.at(i).lastName << std::endl;
        Edge name {s, i+1, instructors.at(i).maxCourses};
        input_vector.push_back(name);
        match.insert(std::pair<std::string, unsigned>(instructors.at(i).lastName, i+1));
        match_reversed.insert(std::pair<unsigned, std::string>(i+1, instructors.at(i).lastName));
        // for (auto x : input_vector) 
        // {
        //     std::cout << x.from << ' ' << x.to << ' ' << x.weight << std::endl;
        // }
    }
    for (unsigned i = 0; i < courses.size(); i++)
    {
        //std::cout << i+instructors.size()+1 << " : " << courses.at(i) << std::endl;
        Edge name {i+(unsigned)instructors.size()+1, t, 1};
        input_vector.push_back(name);
        match.insert(std::pair<std::string, unsigned>(courses.at(i), i+instructors.size()+1));
        match_reversed.insert(std::pair<unsigned, std::string>(i+instructors.size()+1, courses.at(i)));
        // for (auto x : input_vector) 
        // {
        //     std::cout << x.from << ' ' << x.to << ' ' << x.weight << std::endl;
        // }
    }
    for (unsigned i = 0; i < instructors.size(); i++)
    {
        for (unsigned j = 0; j < instructors.at(i).preferences.size(); j++)
        {
            Edge name {match[instructors.at(i).lastName], match[instructors.at(i).preferences.at(j)], 1};
            input_vector.push_back(name);
        }
    }

    std::vector<Edge> course_assigned = solveNetworkFlow(input_vector, t+1);
    // std::cout << "------map-------" << std::endl;
    // std::map<std::string, unsigned>::iterator it = match.begin();
    // while (it != match.end()) {
    //     std::cout << it->first << ' ' << it->second << std::endl;
    //     it ++;
    // }

    // std::cout << "------result-------" << std::endl;
    // for (auto x : course_assigned) 
    //     {
    //         std::cout << x.from << ' ' << x.to << ' ' << x.weight << std::endl;
    //     }

    for (unsigned i = 0; i < course_assigned.size(); i++)
    {
        if (course_assigned.at(i).from != s && course_assigned.at(i).to != t)
        {
            unsigned prof = course_assigned.at(i).from -1;
            std::string course = match_reversed[course_assigned.at(i).to];
            instructors[prof].assignedCourses.push_back(course);
        }
    }
    
}
