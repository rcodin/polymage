#include <Python.h>
#include <iostream>
#include <unordered_map> 
#include <unordered_set>
#include <set>
#include <algorithm>
#include <queue>
#include <regex>
#include <fstream>
#include <string>
#include <boost/multiprecision/cpp_int.hpp> 
#include <functional> 
#include <cctype>
#include <locale>
#include <math.h>
#include <map>
#include <numeric>

#define PRINT_DEBUG (x) if (DEBUG) { x; }
bool const DEBUG = true;
bool const checkForAssertions = false;
bool const checkPythonExceptions = true;

#define __MCASTLE1_SERVER__

#ifdef __PERSONAL_COMPUTER__
static int const IMAGE_ELEMENT_SIZE = 4; //Each image element is of 4 bytes
static int const L2_CACHE_SIZE = 256*1024;  //The L2_CACHE_SIZE is 256KB
static int const N_CORES = 4; //Number of Cores is 4
#endif

#ifdef __MCASTLE1_SERVER__
static int const IMAGE_ELEMENT_SIZE = 4; //Each image element is of 4 bytes
static int const L2_CACHE_SIZE = 512*1024;  //The L2_CACHE_SIZE is 512KB
static int const N_CORES = 16; //Number of Cores is 16
static int const L1_CACHE_SIZE = 32*1024;
#endif

#ifdef __MCASTLE2_SERVER__
static int const IMAGE_ELEMENT_SIZE = 4; //Each image element is of 4 bytes
static int const L2_CACHE_SIZE = 2*512*1024;  //The L2_CACHE_SIZE is 512KB
static int const N_CORES = 16; //Number of Cores is 16
#endif

//To Compile g++ -shared -I/usr/include/python3.4m/ -lpython3 -o optgrouping.so optgrouping_incremental.cpp -fPIC -O3 -std=c++11 -lboost_system
using boost::multiprecision::uint128_t;
using namespace boost;
using namespace std;
uint64_t logMaxChildren = 3;
const uint64_t log_increment_factor = 1;
PyObject* reduction_cls;
PyObject* small_comps;
PyObject* comp_size_map;
PyObject* tstencil_cls;
PyObject* pipeline;
PyObject* pygroup_topological_order;
PyObject* pygroup_dim_reuse;
PyObject* pylive_size;
PyObject* pydim_size;
PyObject* storage_mapping_get_dim_size;
PyObject* cls_Storage;
PyObject* get_overlapping_size_func;

long grp_size;
uint64_t globalID = -1, runningTime = 0, maxID = 0, foundInDict = 0;

struct uint128Hasher
{
    uint64_t operator()(const uint128_t i) const
    {
        //printf ("hash_id %ld\n", g->hashID());
        return (uint64_t)i;
    }
};

struct uint128Comparer
{
    uint64_t operator()(const uint128_t i, const uint128_t j) const
    {
        //printf ("hash_id %ld\n", g->hashID());
        return i < j;
    }
};

struct PyObjectHasher
{
    uint64_t operator() (const PyObject* obj) const
    {
        return (uint64_t)obj;
    }
};

struct PyObjectPred
{
    bool operator() (PyObject* obj1, PyObject* obj2) const
    {
        return (bool)PyObject_RichCompareBool (obj1, obj2, Py_EQ);
    }
};

class OptGroup;

class Group 
{
public:
    struct GroupComparer
    {
        bool operator()(const Group* g1, const Group* g2) const
        {
            return g1->hashID() < g2->hashID();
        }
    };
    
    struct GroupHasher
    {
        uint64_t operator()(const Group* g1) const
        {
            return (uint64_t)g1->hashID();
        }
    };
    
protected:
    uint64_t id;
    uint128_t hash_id;
    vector<uint128_t> next_groups_int;
    vector<uint128_t> prev_groups_int;
    bool is_last;
    PyObject* pyGroup;
    PyObject* comps;
    
    set<Group*, Group::GroupComparer> next_groups, prev_groups;
    set<OptGroup*, Group::GroupComparer> parent_groups;
    
public:
    static std::unordered_map<uint128_t, Group*, uint128Hasher> hashIDToGroup;
    
    Group ()
    {
        setHashID ();
        pyGroup = NULL;
        comps = NULL;
    }
    
    Group (PyObject* pyGroup, PyObject* comps)
    {
        this->pyGroup = pyGroup;
        this->comps = comps;
        setHashID ();
    }
    
    Group (uint128_t hash_id)
    {
        setHashID (hash_id);
        pyGroup = NULL;
        comps = NULL;
    }
    
    Group (set<Group*, GroupComparer>& next_groups, set<Group*, GroupComparer>& prev_groups)
    {
        this->next_groups = next_groups;
        this->prev_groups = prev_groups;
        if (nextGroups().size() == 0)
            is_last = true;
        else
            is_last = false;
        
        id = ++globalID;
        setHashID ();
    }
    
    void setHashID ()
    {
        id = ++globalID;
        if (id < maxID)
        {
            hash_id = 1L;
            hash_id = (hash_id << id);
            hashIDToGroup[hash_id] = this;
        }
    }
    
    PyObject* getPyGroup ()
    {
        return pyGroup;
    }
    
    PyObject* getComps ()
    {
        return comps;
    }
    
    PyObject* getCompAtIndex (int i)
    {
        return PyList_GetItem (comps, i);
    }
    
    set<OptGroup*, GroupComparer>& parentGroups()
    {
        return parent_groups;
    }
    
    set<Group*, GroupComparer>& nextGroups()
    {
        return next_groups;
    }
    
    set<Group*, GroupComparer>& prevGroups ()
    {
        return prev_groups;
    }
    
    void setHashID (uint128_t _hash_id)
    {
        hash_id = _hash_id;
    }
    
    uint128_t hashID () const
    {
        return hash_id;
    }
    
    bool isSingleGroup () const
    {
        if (hash_id && (!(hash_id & (hash_id-1))))
            return true;
        
        return false;
    }
    
    vector<uint128_t>* prevGroupsHashID ()
    {
        if (prev_groups_int.size () != 0)
        {
            return &prev_groups_int;
        }
        
        for (auto it = prevGroups().begin(); it != prevGroups().end(); ++it)
        {
            prev_groups_int.push_back ((*it)->hashID());
        }
         
        return &prev_groups_int;
    }
    
    vector<uint128_t>* nextGroupsHashID ()
    {
        if (next_groups_int.size () != 0)
        {
            return &next_groups_int;
        }
        
        for (auto it = nextGroups().begin(); it != nextGroups().end(); ++it)
        {
            next_groups_int.push_back ((*it)->hashID());
        }
         
        return &next_groups_int;
    }
    
    
    bool operator==(const Group& other) const
    {
        return other.hash_id == hash_id;
    }
};

class OptGroup : public Group
{
protected:
    set<OptGroup*, Group::GroupComparer> children;
    
    OptGroup (uint128_t hash_id) : Group (hash_id)
    {
        //setChildren ();
    }

public:   
    static std::unordered_map<uint128_t, OptGroup*, uint128Hasher> optHashIDToGroup;
    static std::unordered_map<uint128_t, OptGroup*, uint128Hasher> parentOptGroupsToHashID;
    static vector <OptGroup*> vectorOptGroups;
    
    //Find if removing a child will make graph disconnected
    //If true, then fill the vector components with 2 disconnected components
    bool isDisconnected2 (OptGroup* toRemove, std::vector< std::vector <OptGroup*> >& components)
    {
        //Using Kosaraju’s DFS to find if removing the child would make
        //the group disconnected
        
        ////std::cout<<"isDisconnected "<< std::bitset<41>((uint64_t)(hashID())) << " toremove " <<std::bitset<41>((uint64_t)toRemove->hashID()) << std::endl;
        stack<OptGroup*> s;
        unordered_map<uint128_t, bool, uint128Hasher> visited;
        OptGroup* start = NULL;
        bool _isDisconnected = false;
        
        for (auto it = children.begin(); it != children.end (); it++)
        {
            if ((*it)->hashID() != toRemove->hashID())
                start = *it;
            
            visited[(*it)->hashID ()] = false;
        }
        
        ////std::cout<<"start " << std::bitset<41>((uint64_t)(start->hashID()))<<std::endl;
        s.push (start);
        
        while (s.size () != 0)
        {
            OptGroup* g = s.top ();
            s.pop ();
            
            visited[g->hashID()] = true;
            
            for (auto it = g->nextGroups ().begin (); 
                 it != g->nextGroups().end (); it++)
            {
                if (((*it)->hashID() & hashID ()) == (*it)->hashID () &&
                    (*it)->hashID () != toRemove->hashID () && 
                    visited[(*it)->hashID()] == false)
                {
                    s.push ((OptGroup*)*it);
                }
            }
        }
        
        unordered_map<uint128_t, bool, uint128Hasher> reverseVisited;
        
        for (auto it = children.begin(); it != children.end (); it++)
        {
            reverseVisited[(*it)->hashID ()] = false;
        }
        
        s.push (start);
        
        while (s.size () != 0)
        {
            OptGroup* g = s.top ();
            s.pop ();
            
            reverseVisited[g->hashID ()] = true;
            
            for (auto it = g->prevGroups ().begin (); 
                 it != g->prevGroups().end (); it++)
            {
                if (((*it)->hashID() & hashID ()) == (*it)->hashID () &&
                    (*it)->hashID () != toRemove->hashID () &&
                    reverseVisited[(*it)->hashID()] == false)
                {
                    s.push ((OptGroup*)*it);
                }
            }
        }
        
        for (auto it = children.begin (); it != children.end (); it++)
        {
            if ((*it)->hashID() != toRemove->hashID ())
            {            
                if (visited[(*it)->hashID ()] == false &&
                    reverseVisited[(*it)->hashID()] == false)
                {
                    _isDisconnected = true;
                    components[1].push_back (*it);
                }
                else
                {
                    components[0].push_back (*it);
                }
            }   
        }
        
        return _isDisconnected;
    }
    
    bool allPrevGroupsAreNotChild (OptGroup* group)
    {        
        for (auto it = group->prevGroups().begin (); it != group->prevGroups().end(); it++)
        {
            if (((*it)->hashID () & hashID()) == (*it)->hashID ())
            {
                return false;
                break;
            }
        }  
        
        return true;
    }
    
    bool allNextGroupsAreNotChild (OptGroup* group)
    {        
        for (auto it = group->nextGroups().begin (); it != group->nextGroups().end(); it++)
        {
            if (((*it)->hashID () & hashID()) == (*it)->hashID ())
            {
                return false;
                break;
            }
        }  
        
        return true;
    }
    
    bool allPrevGroupsAreChild (OptGroup* group)
    {        
        for (auto it = group->prevGroups().begin (); it != group->prevGroups().end(); it++)
        {
            if (((*it)->hashID () & hashID()) != (*it)->hashID ())
            {
                return false;
                break;
            }
        }  
        
        return true;
    }
    
    bool allNextGroupsAreChild (OptGroup* group)
    {        
        for (auto it = group->nextGroups().begin (); it != group->nextGroups().end(); it++)
        {
            if (((*it)->hashID () & hashID()) != (*it)->hashID ())
            {
                return false;
                break;
            }
        }  
        
        return true;
    }
    
    bool isDisconnected (OptGroup* toRemove, vector<unordered_set <OptGroup*> >& components)
    {
        //If either all of toRemove's previous or next groups are not in this group
        //Then graph will not become disconnected.
        ////std::cout<<"isDisconnected allPrevGroupsNot " << allPrevGroupsAreNotChild (toRemove) << "  " << allNextGroupsAreNotChild (toRemove) << std::endl;
        if (//allPrevGroupsAreNotChild (toRemove) ||
            allNextGroupsAreNotChild (toRemove))
            return false;

        vector<OptGroup*> source_vertices;
        unordered_map<uint128_t, bool, uint128Hasher> visited;
         //vector of vertices reachable from one vertex using nextGroups
        unordered_map<OptGroup*, unordered_set <OptGroup*> > reachable_vertices;
        //vector <unordered_set <OptGroup*> > _components;
        
        for (auto it = children.begin(); it != children.end (); it++)
        {
            stack<OptGroup*> s;
            
            if ((*it)->hashID() == toRemove->hashID())
                continue;
            
            unordered_set <OptGroup*> reachable_set;
            
            s.push ((OptGroup*)*it);
            
            while (s.size () != 0)
            {
                OptGroup* g = s.top ();
                s.pop ();
                
                for (auto it2 = g->nextGroups ().begin (); 
                    it2 != g->nextGroups().end (); it2++)
                {
                    if (((*it2)->hashID() & hashID ()) == (*it2)->hashID () &&
                        (*it2)->hashID () != toRemove->hashID ())
                    {
                        reachable_set.insert ((OptGroup*)*it2);
                        s.push ((OptGroup*)*it2);
                    }
                }
            }
            
            reachable_vertices [*it] = reachable_set;
        }
        
        ////std::cout<< "isDisconnected reachable_vertices: " << std::endl;
        for (auto it = reachable_vertices.begin(); 
             it != reachable_vertices.end(); it++)
        {
            ////std::cout<< std::bitset<41>((uint64_t)it->first->hashID ())<<std::endl;
            
            for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
            {
                ////std::cout<< "    " << std::bitset <41>((uint64_t) (*it2)->hashID ())<< std::endl;
            }
        }
        
        unordered_map <OptGroup*, int> group_to_set_index;
        
        for (auto it = reachable_vertices.begin(); 
             it != reachable_vertices.end(); it++)
        {
            unordered_set <OptGroup*>* foundInSet = NULL;
            int setIndex = 0;
            
            for (auto it2 = components.begin (); it2 != components.end(); it2++)
            {
                if ((*it2).find (it->first) != (*it2).end ())
                {
                    foundInSet = &(*it2);
                    //setIndex = it2 - components.begin();
                    break;
                }
            }
            
            if (foundInSet == NULL)
            {
                for (auto it2 = it->second.begin(); it2 != it->second.end (); it2++)
                {
                    if (group_to_set_index.find ((OptGroup*)*it2) != 
                        group_to_set_index.end ())
                    {
                        foundInSet = &components [group_to_set_index [(OptGroup*)*it2]];
                        break;
                    }
                }
            }
            
            if (foundInSet != NULL)
            {
                for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
                {
                    foundInSet->insert ((OptGroup*)*it2);
                }
                
                foundInSet->insert (it->first);
            }
            else
            {                
                unordered_set<OptGroup*> s = it->second;
                s.insert (it->first);
                components.push_back (s);
                group_to_set_index [it->first] = components.size () - 1;
            }
        }
        
        if (components.size () == 1)
            return false;
            
        return true;
    }
    
    bool createsCycle (OptGroup* toRemove, std::vector< std::vector <OptGroup*> >& components)
    {
        //If toRemove's previous group's (in the parent) can reach toRemove's next group's (in the parent)
        //by not passing through toRemove, then there will be cycle after removing toRemove
        
        vector <OptGroup*> prev_children;
        vector <OptGroup*> next_children;
        bool is_createsCycle = false;
        
        for (auto it = toRemove->prevGroups ().begin ();
             it != toRemove->prevGroups ().begin (); it++)
        {
            if ((hashID()& (*it)->hashID()) == (*it)->hashID())
                prev_children.push_back ((OptGroup*)*it);
        }
        
        for (auto it = toRemove->nextGroups ().begin ();
             it != toRemove->nextGroups ().begin (); it++)
        {
            if ((hashID()& (*it)->hashID()) == (*it)->hashID())
                next_children.push_back ((OptGroup*)*it);
        }
        
        vector<OptGroup*> source_vertices;
        unordered_map<uint128_t, bool, uint128Hasher> visited;
        stack<OptGroup*> s;
        
        for (auto it = prev_children.begin(); it != prev_children.end (); it++)
        {
            s.push (*it);
        }
        
        while (s.size () != 0)
        {
            OptGroup* g = s.top ();
            s.pop ();
            
            if (toRemove->nextGroups ().find (g) != 
                toRemove->nextGroups ().end())
            {
                is_createsCycle = true;
                break;
            }
            
            for (auto it = g->nextGroups ().begin ();
                 it != g->nextGroups().end (); it++)
            {
                if (((*it)->hashID() & hashID ()) == (*it)->hashID () &&
                    (*it)->hashID () != toRemove->hashID ())
                {
                    s.push ((OptGroup*)*it);
                }
            }
        }
        
        ////std::cout<< "is_createsCycle "<<is_createsCycle << std::endl;
        if (!is_createsCycle)
            return false;
        
        while (s.size () != 0)
            s.pop ();
        
        for (auto it = next_children.begin(); it != next_children.end (); it++)
        {
            s.push (*it);
        }
        
        while (s.size () != 0)
        {
            OptGroup* g = s.top ();
            s.pop ();
            
            visited[g->hashID()] = true;
            
            for (auto it = g->nextGroups ().begin (); 
                 it != g->nextGroups().end (); it++)
            {
                if (((*it)->hashID() & hashID ()) == (*it)->hashID () &&
                    (*it)->hashID () != toRemove->hashID ())
                {
                    s.push ((OptGroup*)*it);
                }
            }
        }
        
        for (auto it = children.begin (); it != children.end (); it++)
        {
            if ((*it)->hashID() != toRemove->hashID ())
            {            
                if (visited[(*it)->hashID ()] == false)
                {
                    components[1].push_back (*it);
                }
                else
                {
                    components[0].push_back (*it);
                }
            }   
        }
        
        if (components[1].size () == 0)
            return false;
            
        return true;
    }
    
    bool minParentEdgeToChild (OptGroup* toRemove, std::vector< unordered_set <OptGroup*> >& components, 
                       OptGroup* minParentGroup)
    {
        if (!minParentGroup)
        {
            printf ("minParentEdgeToChild: minParentGroup is NULL\n");
            assert (false);
        }
            
        //If a child of min parent group has an edge to a child of this group
        //then break this group as it can create cycles
        
        bool toReturn = false;
        vector<OptGroup*> childrenToDisconnect;
        
        for (auto it = minParentGroup->childGroups ().begin(); 
             it != minParentGroup->childGroups ().end(); it++)
        {
            ////std::cout<< "minParentGroup child " << std::bitset <41>((uint64_t)(*it)->hashID())<<std::endl;
            for (auto it2 = (*it)->nextGroups ().begin(); 
                 it2 != (*it)->nextGroups ().end(); it2++)
            {
                ////std::cout<< "child next "<< std::bitset <41>((uint64_t)(*it2)->hashID())<<std::endl;
                if (((*it2)->hashID () & hashID ()) == (*it2)->hashID ())
                {
                    ////std::cout<< "is parent's child "<<std::endl;
                    childrenToDisconnect.push_back ((OptGroup*)(*it2));
                    toReturn = true;
                }
            }
        }
        
        ////std::cout<< "minParentEdgeToChild " << toReturn << std::endl;
        
        ////std::cout<< "childrenToDisconnect "<< childrenToDisconnect.size () << std::endl;
        for (auto it = childrenToDisconnect.begin(); it != childrenToDisconnect.end (); it++)
        {
            ////std::cout<< std::bitset<41>((uint64_t)(*it)->hashID ())<< std::endl;
        }
        
        if (!toReturn)
            return toReturn;
        
        vector<OptGroup*> source_vertices;
        unordered_map<uint128_t, bool, uint128Hasher> visited;
        unordered_map<OptGroup*, unordered_set <OptGroup*> > reachable_vertices;
        
        components.clear ();
        
        for (auto it = children.begin(); it != children.end (); it++)
        {
            visited [(*it)->hashID ()] = false;
        }
        
        for (auto it = childrenToDisconnect.begin(); 
             it != childrenToDisconnect.end (); it++)
        {
            stack<OptGroup*> s;
            unordered_set <OptGroup*> reachable_set;
            s.push (*it);
            
            while (s.size () != 0)
            {
                OptGroup* g = s.top ();
                s.pop ();
                visited[g->hashID()] = true;
                
                for (auto it2 = g->nextGroups ().begin ();
                    it2 != g->nextGroups().end (); it2++)
                {
                    if (((*it2)->hashID() & hashID ()) == (*it2)->hashID () &&
                        (*it2)->hashID () != toRemove->hashID ())
                    {
                        reachable_set.insert ((OptGroup*)*it2);
                        s.push ((OptGroup*)*it2);
                    }
                }
            }
            
            reachable_vertices [*it] = reachable_set;
        }

        components.push_back (unordered_set <OptGroup*> ());
        
        for (auto it = children.begin (); it != children.end (); it++)
        {
            ////std::cout<< "isvisited " << std::bitset<41>((uint64_t)(*it)->hashID ()) << " " << visited[(*it)->hashID ()]<<std::endl;
            if (visited[(*it)->hashID ()] == false)
            {            
                components[0].insert (*it);
            }   
        }
        
        if (components [0].size () == 0)
            components.clear ();
        
        ////std::cout<< "minParentEdgeToChild reachable_vertices: " << std::endl;
        for (auto it = reachable_vertices.begin(); 
             it != reachable_vertices.end(); it++)
        {
            ////std::cout<< std::bitset<41>((uint64_t)it->first->hashID ())<<std::endl;
            
            for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
            {
                ////std::cout<< "    " << std::bitset <41>((uint64_t) (*it2)->hashID ())<< std::endl;
            }
        }
        
        ////std::cout<< "components EARLIER 0" << std::endl;
        for (auto it = components.begin(); it != components.end(); it++)
        {
            ////std::cout<< "component " << (it)-components.begin() << std::endl;
            for (auto it2 = (*it).begin (); it2 != (*it).end(); it2++)
            {
                ////std::cout<< "   " << std::bitset <41>((uint64_t)(*it2)->hashID()) << std::endl;
            }
        }
        
        for (auto it = reachable_vertices.begin(); 
             it != reachable_vertices.end(); it++)
        {
            unordered_set <OptGroup*>* foundInSet = NULL;
            
            for (auto it2 = components.begin (); it2 != components.end(); it2++)
            {
                if ((*it2).find (it->first) != (*it2).end ())
                {
                    foundInSet = &(*it2);
                    break;
                }
            }
            
            if (foundInSet != NULL)
            {
                for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
                {
                    foundInSet->insert ((OptGroup*)*it2);
                }
                
                foundInSet->insert (it->first);
            }
            else
            {
                components.push_back (it->second);
                (*(components.end () - 1)).insert (it->first);
            }
        }
        
        ////std::cout<< "components" << std::endl;
        for (auto it = components.begin(); it != components.end(); it++)
        {
            ////std::cout<< "component " << (it)-components.begin() << std::endl;
            for (auto it2 = (*it).begin (); it2 != (*it).end(); it2++)
            {
                ////std::cout<< "   " << std::bitset <41>((uint64_t)(*it2)->hashID()) << std::endl;
            }
        }
        
        if (components.size () == 1)
            return false;
            
        return true;
    }
    
    void removeChild (uint128_t _hash_id, OptGroup* minParentGroup)
    {
        vector <unordered_set <OptGroup*> > components;
        std::vector <OptGroup*> childrenToRemove;
        uint128_t prev_hash_id = hashID();

        ////std::cout<< "removing from parent " << std::bitset<41>((uint64_t)(hashID()))<<std::endl;
        if (this->hashID() == 0 || _hash_id == 0)
            return;
        
        OptGroup* toRemove = OptGroup::optHashIDToGroup[_hash_id];
        childrenToRemove.push_back (toRemove);
        ////std::cout<<"toRemove is " << std::bitset<41>((uint64_t)toRemove->hashID()) << " parent group address " <<(uint64_t)this<<std::endl;
        //TODO: create functions for setting/getting values in parentOptGroupToHashID map and optGroupToHashID map
        //Also add assertions in it that the group being set/get is equal to the hash id key.
        for (auto it = prevGroups().begin(); it != prevGroups().end();
             it++)
        {
            ////std::cout<<"before removing prevgroup " << std::bitset<41>((uint64_t)(*it)->hashID())<< "  " << (uint64_t)*it <<std::endl;
        }
        
        for (auto it = nextGroups().begin(); it != nextGroups().end();
             it++)
        {
            ////std::cout<<"before removing nextGroups " << std::bitset<41>((uint64_t)(*it)->hashID()) << "  " << (uint64_t)*it<<std::endl;
        }
        
        if (_hash_id != hashID () &&
            (isDisconnected (childrenToRemove[0], 
            components) || minParentEdgeToChild (childrenToRemove[0], 
            components, minParentGroup)))
        {
            if (components.size() > 1)
            {
                ////std::cout<<"isDisconnected is true " << components.size () << std::endl;
                auto components_it = components.begin ();
                
                components_it++; //the first one is current group
                while (components_it != components.end ())
                {
                    unordered_set <OptGroup*>& connected_component = *components_it;
                    ////std::cout<< "connected_component size " << connected_component.size () << std::endl;
                    uint128_t newHashID = (*connected_component.begin())->hashID();
                    
                    for (auto it = connected_component.begin (); it != connected_component.end ();
                        it++)
                    {
                        newHashID = newHashID | (*it)->hashID ();
                        ////std::cout<<"removing "<< std::bitset<41>((uint64_t)(*it)->hashID ()) << std::endl;
                        children.erase (OptGroup::optHashIDToGroup[(*it)->hashID ()]);
                        childrenToRemove.push_back (*it);
                    }
                
                    ////std::cout<< "newHashID " << std::bitset<41>((uint64_t)newHashID) << " in parentopt " << inParentOptGroupsToHashID (newHashID) << std::endl;
                    OptGroup* newGroup;
                    
                    if (!inParentOptGroupsToHashID (newHashID))
                    {
                        newGroup = OptGroup::createParentOptGroup (newHashID);
                        vector <Group*> toRemove;
                    
                        //Find next and prev groups for new group
                        
                        for (auto it2 = connected_component.begin (); 
                             it2 != connected_component.end (); it2++)
                        {
                            for (auto it3 = (*it2)->nextGroups ().begin (); 
                                 it3 != (*it2)->nextGroups ().end (); it3++)
                            {
                                for (auto it4 = (*it3)->parentGroups ().begin ();
                                     it4 != (*it3)->parentGroups ().end ();
                                     it4++)
                                {
                                    if ((*it4)->hashID () != newGroup->hashID ())
                                    {
                                        ////std::cout<<"newgroup nextGroup " << std::bitset<41>((uint64_t)(*it4)->hashID())<<std::endl;
                                        newGroup->nextGroups ().insert (*it4);
                                        (*it4)->prevGroups (). insert (newGroup);
                                    }
                                }
                            }
                        }
                        
                        for (auto it2 = connected_component.begin (); 
                             it2 != connected_component.end (); it2++)
                        {
                            for (auto it3 = (*it2)->prevGroups ().begin (); 
                                 it3 != (*it2)->prevGroups ().end (); it3++)
                            {
                                for (auto it4 = (*it3)->parentGroups ().begin ();
                                     it4 != (*it3)->parentGroups ().end ();
                                     it4++)
                                {
                                    if ((*it4)->hashID () != newGroup->hashID ())
                                    {
                                        ////std::cout<<"newgroup prevGroups " << std::bitset<41>((uint64_t)(*it4)->hashID())<<std::endl;
                                        newGroup->prevGroups ().insert (*it4);
                                        (*it4)->nextGroups (). insert (newGroup);
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        newGroup = parentOptGroupsToHashID[newHashID];
                        
                       // //std::cout<< "newGroup in parent " << newGroup <<  std::endl;
                    }
                    
                    //nextGroups ().insert (newGroup);
                    //newGroup->prevGroups (). insert (this);
                    
                    components_it++;
                }
            }
        }
        
        /*if (hash_id == 0)
        {
            //OptGroup is removed
            for (auto it = children.begin(); it != children.end(); it++)
            {
                (*it)->parentGroups().erase (this);
            }
        }
        else*/
        
        for (auto const &childToRemove : childrenToRemove)
        {
            vector<OptGroup*> prevGroupsToRemove, nextGroupsToRemove;
            
            ////std::cout<< std::bitset<41>((uint64_t)(hash_id))<< "  " << std::bitset<41>((uint64_t)(childToRemove)->hashID())<<std::endl;
            ////std::cout<< "contains parent " << ((childToRemove)->parentGroups().find(this) != (childToRemove)->parentGroups().end()) << std::endl;
            (childToRemove)->parentGroups().erase (this);
            
            for (auto _it = (childToRemove)->parentGroups().begin();
                _it != (childToRemove)->parentGroups().end ();
                _it++)
            {
                //std::cout<<"after removing parent group:" << std::bitset<41>((uint64_t)(*_it)->hashID ()) << std::endl;
            }
            
            children.erase (childToRemove);
            
            //Remove those links from prevGroups and nextGroups to this group
            //which are due to children in childrenToRemove
            
            //std::cout<<"remove prev next for " << std::bitset<41>((uint64_t)(childToRemove)->hashID()) << std::endl;
            for (auto it2 = prevGroups().begin(); it2 != prevGroups().end(); it2++)
            {
                OptGroup* prevGrp = (OptGroup*)(*it2);
                
                for (auto it3 = childToRemove->prevGroups().begin (); 
                     it3 != childToRemove->prevGroups().end(); it3++)
                {
                    if ((*it3)->parentGroups().find (prevGrp) != (*it3)->parentGroups().end())
                    {
                        bool canRemove = true;
                        for (auto it4 = children.begin(); it4 != children.end();
                            it4++)
                        {
                            for (auto it5 = (*it4)->prevGroups().begin(); 
                                 it5 != (*it4)->prevGroups().end(); it5++)
                            {
                                //Remove parent group only if it is the parent of only
                                //prev group of only one child
                                if ((*it5)->parentGroups().find (prevGrp) != (*it5)->parentGroups().end())
                                {
                                    canRemove = false;
                                    break;
                                }
                            }
                        }
                        
                        if (canRemove)
                            prevGroupsToRemove.push_back (prevGrp);
                    }
                }
            }
            
            for (auto it2 = nextGroups().begin(); it2 != nextGroups().end(); it2++)
            {
                OptGroup* nextGrp = (OptGroup*)(*it2);
                
                for (auto it3 = childToRemove->nextGroups().begin (); 
                     it3 != childToRemove->nextGroups().end(); it3++)
                {
                    if ((*it3)->parentGroups().find (nextGrp) != 
                        (*it3)->parentGroups().end())
                    {
                        bool canRemove = true;
                        for (auto it4 = children.begin(); it4 != children.end();
                            it4++)
                        {
                            for (auto it5 = (*it4)->nextGroups().begin(); 
                                 it5 != (*it4)->nextGroups().end(); it5++)
                            {
                                //Remove parent group only if it is the parent of only
                                //prev group of only one child
                                if ((*it5)->parentGroups().find (nextGrp) != (*it5)->parentGroups().end())
                                {
                                    canRemove = false;
                                    break;
                                }
                            }
                        }
                        
                        if (canRemove)
                            nextGroupsToRemove.push_back (nextGrp);
                    }
                }
            }
            
            for (auto it2 = prevGroupsToRemove.begin();
                 it2 != prevGroupsToRemove.end(); it2++)
            {
                //std::cout<<"remove prev " << std::bitset<41>((uint64_t)(*it2)->hashID()) << std::endl;
                prevGroups().erase (*it2);
                (*it2)->nextGroups().erase (this);
            }
            
            for (auto it2 = nextGroupsToRemove.begin();
                 it2 != nextGroupsToRemove.end(); it2++)
            {
                //std::cout<<"remove next " << std::bitset<41>((uint64_t)(*it2)->hashID()) << std::endl;
                nextGroups().erase (*it2);
                (*it2)->prevGroups().erase (this);
            }
        }
        
        //std::cout<< "number of children " << children.size () << std::endl;
        //Clear prevGroups and nextGroups
        
        prevGroups().clear ();
        nextGroups().clear ();
        
        //Update prevGroup and nextGroup with updated children set
        for (auto const &child : children)
        {
            //Add parentGroups of prevGroups of child as prevGroups
            for (auto const &child_prev : child->prevGroups ())
            {
                if (children.find ((OptGroup*)child_prev) == children.end ())
                    prevGroups().insert (child_prev->parentGroups().begin(),
                                         child_prev->parentGroups().end());
            }
            
            //Add parentGroups of nextGroups of child as nextGroups
            for (auto const &child_next : child->nextGroups ())
            {
                if (children.find ((OptGroup*)child_next) == children.end ())
                    nextGroups().insert (child_next->parentGroups().begin(),
                                         child_next->parentGroups().end());
            }
        }
        
        for (auto it = prevGroups().begin(); it != prevGroups().end();
             it++)
        {
            //std::cout<<"after removing prevgroup " << std::bitset<41>((uint64_t)(*it)->hashID())<<std::endl;
        }
        
        for (auto it = nextGroups().begin(); it != nextGroups().end();
             it++)
        {
            //std::cout<<"after removing nextGroups " << std::bitset<41>((uint64_t)(*it)->hashID())<<std::endl;
            for (auto it2 = (*it)->nextGroups().begin(); it2 != (*it)->nextGroups().end();
             it2++)
            {
                //std::cout<<"next " << std::bitset<41>((uint64_t)(*it2)->hashID())<<std::endl;
            }
        }
        
        //Update all sets where this group is contained because OptGroup
        //is hashed on the base of hash id. Changing hash id of group without 
        //updating there could (or will?) produce error
        for (auto it = children.begin(); it != children.end(); it++)
        {
            //Remove from parentGroup of children
            (*it)->parentGroups().erase (this);
        }
        
        for (auto it = nextGroups().begin(); it != nextGroups().end(); it++)
        {
            //Remove from prevGroup of next groups
            (*it)->prevGroups().erase (this);
        }
        
        for (auto it = prevGroups().begin(); it != prevGroups().end(); it++)
        {
            //Remove from nextGroup of prev groups
            (*it)->nextGroups().erase (this);
        }
        
        //Update Hash ID now because set of parents use hash id as its hash id
        for (auto it = childrenToRemove.begin (); it != childrenToRemove.end ();
             it++)
        {
            hash_id &= ~(*it)->hashID ();
        }
        
        if (hash_id != 0)
        {
            for (auto it = children.begin(); it != children.end(); it++)
            {
                //Insert to parentGroup of children
                (*it)->parentGroups().insert (this);
            }
            
            for (auto it = nextGroups().begin(); it != nextGroups().end(); it++)
            {
                //Insert to prevGroup of next groups
                (*it)->prevGroups().insert (this);
            }
            
            for (auto it = prevGroups().begin(); it != prevGroups().end(); it++)
            {
                //Insert to nextGroup of prev groups
                (*it)->nextGroups().insert (this);
            }
        }
        
        for (auto it = nextGroups().begin(); it != nextGroups().end();
             it++)
        {
            //std::cout<<"after removing121212 nextGroups " << std::bitset<41>((uint64_t)(*it)->hashID())<<std::endl;
            for (auto it2 = (*it)->nextGroups().begin(); it2 != (*it)->nextGroups().end();
             it2++)
            {
                //std::cout<<"next " << std::bitset<41>((uint64_t)(*it2)->hashID())<<std::endl;
            }
        }
        
        //std::cout<< "new hashid for group " << std::bitset<41>((uint64_t)hashID())<<std::endl;
        if (hash_id != 0)
            parentOptGroupsToHashID [hash_id] = this;
        //std::cout<<"qqqqqqqqqqqqqqqqqq " << std::endl;    
        parentOptGroupsToHashID.erase (prev_hash_id);
        //std::cout<<"qwqwqwq " << std::endl;
        //delete &components[0];
        //delete &components[1];
    }
    
    void setChildrenFromHashID ()
    {        
        for (auto it = optHashIDToGroup.begin (); it != optHashIDToGroup.end (); ++it)
        {
            if (it->second != NULL && it->second != this && (hash_id & it->first) == it->first)
            {
                children.insert (it->second);
                it->second->parentGroups ().insert (this);
            }
        }
    }
    
    set<OptGroup*, GroupComparer>& childGroups ()
    {
        return children;
    }
    
    /*void remove_from_parent ()
    {
        uint64_t bit = 0;
        uint128_t _hash_id = hash_id;
        
        while (_hash_id != 0)
        {
            if (_hash_id & 1L == 1)
            {   
                uint128_t l = 1;
                Group* a = Group::hashIDToGroup[l<<bit];
                for (auto it = a->parentGroups().begin (); it != a->parentGroups().end (); it++)
                {
                    if (*it == this)
                    {
                        a->parentGroups().erase (it);
                        break;
                    }
                }
            }
            
            _hash_id = _hash_id >> 1;
            bit++;
        }
    }*/
    
    inline bool operator==(const OptGroup& other) const
    {
        return other.hash_id == hash_id;
    }
    
    static OptGroup* createOptGroup (uint128_t _hash_id)
    {
        OptGroup* _opt_group = optHashIDToGroup[_hash_id];
        
        if (_opt_group == NULL)
        {
            _opt_group = new OptGroup (_hash_id);
            optHashIDToGroup[_hash_id] = _opt_group;
            vectorOptGroups.push_back (_opt_group);
        }
        
        return (OptGroup*)_opt_group;
    }
    
    static bool inParentOptGroupsToHashID (uint128_t _hash_id)
    {
        return parentOptGroupsToHashID.find(_hash_id) != 
            parentOptGroupsToHashID.end() && parentOptGroupsToHashID [_hash_id] != NULL;
    }
    
    static OptGroup* createParentOptGroup (uint128_t _hash_id)
    {
        OptGroup* _opt_group = parentOptGroupsToHashID[_hash_id];
        
        if (_opt_group == NULL)
        {
            _opt_group = new OptGroup (_hash_id);
            parentOptGroupsToHashID[_hash_id] = _opt_group;
            _opt_group->setChildrenFromHashID ();
        }
        
        return _opt_group;
    }
    
    static vector<uint128_t>* nextGroupsHashID (uint128_t hash_id)
    {
        Group* g;
        
        g = optHashIDToGroup[hash_id];

        if (g != NULL)
            return new vector<uint128_t> (g->nextGroupsHashID()->begin(), 
                                          g->nextGroupsHashID()->end());
                        
        uint64_t bit = 0;
        vector<uint128_t>* next_groups;
        next_groups = new vector<uint128_t> ();
        set <uint128_t, uint128Comparer> _set;
         
        for (auto it = vectorOptGroups.begin (); it != vectorOptGroups.end (); ++it)
        {
            if ((hash_id & (*it)->hashID()) == (*it)->hashID())
            {
                vector<uint128_t>* a = ((Group*)(*it))->nextGroupsHashID();

                for (auto it2 = a->begin(); it2 != a->end(); ++it2)
                {
                    if (((*it2) & hash_id) != (*it2))
                        next_groups->push_back (*it2);
                        //_set.insert (*it2);
                }
            }
        }
        
        /*for (auto it = optHashIDToGroup.begin (); it != optHashIDToGroup.end (); ++it)
        {
            if (it->second != NULL && (hash_id & it->first) == it->first)
            {
                vector<uint128_t>* a = ((Group*)(it->second))->nextGroupsHashID();

                for (auto it2 = a->begin(); it2 != a->end(); ++it2)
                {
                    if (((*it2) & hash_id) != (*it2))
                        next_groups->push_back (*it2);
                        //_set.insert (*it2);
                }
            }
        }*/
        
        for (auto it = _set.begin(); it != _set.end(); it++)
            next_groups->push_back (*it);
          
        return next_groups;
    }
    
    static vector<uint128_t>* prevGroupsHashID (uint128_t hash_id)
    {
        Group* g;
        
        g = optHashIDToGroup[hash_id];

        if (g != NULL)
            return new vector<uint128_t> (g->prevGroupsHashID()->begin(), 
                                          g->prevGroupsHashID()->end());
                        
        uint64_t bit = 0;
        vector<uint128_t>* prev_groups;
        prev_groups = new vector<uint128_t> ();
        
        for (auto it = optHashIDToGroup.begin (); it != optHashIDToGroup.end (); ++it)
        {
            if (it->second != NULL && (hash_id & it->first) == it->first)
            {
                vector<uint128_t>* a = ((Group*)it->second)->prevGroupsHashID();

                for (auto it2 = a->begin(); it2 != a->end(); ++it2)
                {
                    if (((*it2) & hash_id) != (*it2))
                        prev_groups->push_back (*it2);
                }
            }
        }
          
        return prev_groups;
    }
    
    static bool isLiveoutForGroup (uint128_t group, OptGroup* child)
    {        
        for (auto const& next_group: child->nextGroups ())
        {
            if ((next_group->hashID () & group) != next_group->hashID ())
            {
                return true;
            }
        }
        
        if (child->nextGroups ().size () == 0)
            return true;
        
        return false;
    }
    
    static void liveinsForChildInGroup (uint128_t group, OptGroup* child, 
                                        unordered_set<Group*, Group::GroupHasher>& liveins)
    {        
        for (auto const& prev_group: child->prevGroups ())
        {
            if ((prev_group->hashID () & group) != prev_group->hashID ())
            {
                liveins.insert (prev_group);
            }
        }
    }
};

bool optGroupVecFunc (OptGroup* g1, OptGroup* g2)
{
    return (g1->hashID() < g2->hashID());
}


inline uint64_t numberOfOnes (uint128_t hash_id)
{
    uint64_t bit = 0;
    uint64_t nOnes = 0;
        
    while (hash_id != 0)
    {
        if (hash_id & 1L == 1)
        {   
            nOnes ++;
        }
        
        hash_id = hash_id >> 1;
        bit++;
    }
    
    return nOnes;
}

inline uint64_t first_one (uint128_t hash_id)
{
    uint64_t bit = 0;
        
    while (hash_id != 0)
    {
        if (hash_id & 1L == 1)
        {   
            return bit;
        }
        
        hash_id = hash_id >> 1;
        bit++;
    }
    
    return -1;
}

std::unordered_map<uint128_t, uint64_t, uint128Hasher> T;
std::unordered_map<uint128_t, Group*, uint128Hasher> Group::hashIDToGroup;
std::unordered_map<uint128_t, OptGroup*, uint128Hasher> OptGroup::optHashIDToGroup;
std::unordered_map<uint128_t, uint128_t, uint128Hasher> hashIDToOptHashID;
std::unordered_map<PyObject*, Group*> pyToGroup;
std::unordered_map<Group*, PyObject*> groupToPyGroup;
std::unordered_set<PyObject*, PyObjectHasher, PyObjectPred> small_comps_set;
std::unordered_map<uint128_t, OptGroup*, uint128Hasher> OptGroup::parentOptGroupsToHashID;
std::unordered_map<uint128_t, std::vector<uint64_t>, uint128Hasher> optHashIDToTileSizes;
std::unordered_map<PyObject*, PyObject*> pyGroupForComp;
std::unordered_map<PyObject*, PyObject*> comp_to_orig_storage;
std::unordered_map<PyObject*, string> stg_class_to_key;
std::unordered_map<PyObject*, PyObject*> comp_to_func;
std::unordered_map<PyObject*, PyObject*> func_to_typ;
std::unordered_map<PyObject*, long> func_to_ndims;
std::unordered_map<PyObject*, vector<long> > comp_to_max_offsets;
std::unordered_map<PyObject*, vector<PyObject*> > helper_storage_to_dim_storage;
vector<OptGroup*> OptGroup::vectorOptGroups;

uint64_t get_dim_size_for_pygroup (PyObject* group, int dim)
{
    PyObject* dims_size;
    
    dims_size = PyDict_GetItem (pydim_size, group);
    
    return (uint64_t)PyLong_AsLong (PyList_GetItem (dims_size, dim));
}

int get_dim_reuse_for_pygroup (PyObject* group, int dim)
{
    PyObject* dim_reuse;
    
    dim_reuse = PyDict_GetItem (pygroup_dim_reuse, group);
    
    return (int)PyLong_AsLong (PyList_GetItem (dim_reuse, dim));
}

void check_and_print_exception ()
{
    ////std::cout<< "Checking for Exception " << std::endl;
    if (PyErr_Occurred () != NULL)
    {
        //std::cout<<"Exception occurred "<<std::endl;
        PyErr_Print();
        abort ();
    }
}

int get_n_dimensions (PyObject* group)
{
    PyObject* dim_reuse;
    
    dim_reuse = PyDict_GetItem (pygroup_dim_reuse, group);
    //check_and_print_exception ();
    return (int)PyList_Size (dim_reuse);
}

uint64_t get_live_size (PyObject* group)
{
    PyObject* size;
    
    size = PyDict_GetItem (pylive_size, group);
    
    return (uint64_t)PyLong_AsLong (size);
}

int get_max_dimensions_for_group (uint128_t hash_id)
{
    uint64_t bit = 0;
    int max_dim = 0;

    while (hash_id != 0)
    {
        if (hash_id & 1L == 1)
        {
            uint128_t l = 1;
            ////std::cout<< "loop "<< std::bitset<41>((uint64_t)(l << bit)) << std::endl;
            PyObject* pyGroup = Group::hashIDToGroup[l<<bit]->getPyGroup ();
            ////std::cout<< "loop2 "<< std::bitset<41>((uint64_t)(l << bit)) << " " << pyGroup << std::endl;
    
            if (pyGroup != NULL)
            {
                int dim = get_n_dimensions (pyGroup);
                
                if (max_dim < dim)
                    max_dim = dim;
            }
        }
        
        hash_id = hash_id >> 1;
        bit++;
    }
    
    return max_dim;
}

PyObject* getCompForPyGroup (PyObject* pyGroup, int index)
{
    static PyObject* str_comps = Py_BuildValue ("s", "comps");
    PyObject* comps = PyObject_GetAttr (pyGroup, str_comps);
    return PyList_GetItem (comps, 0);
}

char* getPyCompFuncName (PyObject* comp)
{
    static PyObject* str_func = Py_BuildValue ("s", "func");
    static PyObject* str_name = Py_BuildValue ("s", "name");
    PyObject* func = PyObject_GetAttr (comp, str_func);
    PyObject* name = PyObject_GetAttr (func, str_name);
    PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
    return PyBytes_AS_STRING (pyStr);
}

char* getPyGroupName (PyObject* pygroup)
{
    PyObject* name = PyObject_Str (pygroup);
    PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
    return PyBytes_AS_STRING (pyStr);
}

uint64_t getNumberOfInputs (uint128_t hash_id)
{
    uint64_t bit = 0;
    uint64_t nOnes = 0;
    uint64_t input_size = 0;
    uint128_t _hash_id = hash_id;
    
    while (_hash_id != 0)
    {
        if (_hash_id & 1L == 1)
        {
            uint128_t l;
            Group* group;
            PyObject* pyGroup;
            
            l = 1;
            group = Group::hashIDToGroup[l<<bit];
            pyGroup = group->getPyGroup ();
            nOnes++;
            
            if (group->prevGroups().size () == 0)
            {
                input_size++;
            }
            
            for (auto it = group->prevGroups().begin(); 
                 it != group->prevGroups().end (); it++)
            {
                if ((hash_id & (*it)->hashID()) != (*it)->hashID())
                {
                    input_size++;
                }
            }
            
            /*if (false && pyGroup != NULL)
            {
                for (int i = 0; i < get_n_dimensions (pyGroup); i++)
                {
                    //TODO: Optimize it because for every get_dim_reuse_for_pygroup call
                    //dict and list are accessed again. Make them access only once 
                    //and then iterate over list
                    PyObject* comp = PyList_GetItem (comps, 0);
                    PyObject* func = PyObject_GetAttr (comp, str_func);
                    ////std::cout<< "PyGroup " << getPyGroupName (func) << " dim "<< i << " reuse " << get_dim_reuse_for_pygroup (pyGroup, i) << " size " << get_dim_size_for_pygroup (pyGroup, i) << std::endl;
                    dim_reuse [i] += get_dim_reuse_for_pygroup (pyGroup, i)*IMAGE_ELEMENT_SIZE;
                    dim_size [i] += get_dim_size_for_pygroup (pyGroup, i)*IMAGE_ELEMENT_SIZE;
                }
                
                if (max_live_data_size < get_live_size (pyGroup))
                {
                    max_live_data_size = get_live_size (pyGroup);
                }
            }*/
        }
        
        _hash_id = _hash_id >> 1;
        bit++;
    }
    
    return input_size;
}

void get_level_order (uint128_t hash_id, vector <PyObject*>& objs, 
                      const vector <Group*>& vec_groups, 
                      map <Group*, int, Group::GroupComparer>& order)
{
    static PyObject* str_func = Py_BuildValue ("s", "func");
    static PyObject* str_comps = Py_BuildValue ("s", "comps");

    bool change = true;
    
    for (auto it = vec_groups.begin(); it != vec_groups.end (); it ++)
    {
        order [*it] = 0;
    }
    //std::cout<< "34234234 " << std::endl;
    while (change)
    {
        change = false;
        //std::cout<< "change " << change << std::endl;
        for (auto obj = vec_groups.begin (); obj != vec_groups.end (); obj++)
        {
            for (auto prev = (*obj)->prevGroups().begin(); 
                 prev != (*obj)->prevGroups().end(); prev++)
            {
                if (((*prev)->hashID () & hash_id) == (*prev)->hashID ())
                {                
                    if (order.find ((*prev)) != order.end () &&
                        order [(*prev)] >= order [*obj])
                    {
                        order [*obj] = order [(*prev)] + 1;
                        change = true;
                    }
                }
            }
        }
    }
    
    if (DEBUG)
    {
        std::cout<< "optorder is"<<std::endl;
        for (auto it = order.begin(); it != order.end(); it++)
        {
            PyObject* func = PyObject_GetAttr (it->first->getCompAtIndex(0), str_func);
            
            PyObject* name = PyObject_GetAttr (func,
                                            Py_BuildValue ("s", "name"));
            PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
            const char* ss = PyBytes_AS_STRING (pyStr);
                std::cout<< "     "<< ss<< it->second << std::endl;
        }
    }
}

void naive_sched_objs (const map <Group*, int, Group::GroupComparer>& order, unordered_map <Group*, int>& naive_order)
{
    static PyObject* str_func = Py_BuildValue ("s", "func");
    unordered_map <int, vector<Group*>> reverse_map;
    int max_level = 0;
    
    for (auto it = order.begin(); it != order.end(); it++)
    {
        int l = it->second;
        if (reverse_map.find (l) == reverse_map.end())
        {
            reverse_map [l] = std::vector <Group*> ();
            reverse_map [l].push_back (it->first);
        }
        else
        {
            reverse_map [l].push_back (it->first);
        }
        
        max_level = max (max_level, l);
    }
    
    std::cout << "reverse_map: " << std::endl;
    for (auto it = reverse_map.begin (); it != reverse_map.end (); it++)
    {
        std::cout << it->first << std::endl;
        for (auto& it2 : it->second)
        {
            std::cout << "    " <<getPyGroupName (it2->getPyGroup ()) << std::endl;
        }
    }
    
    int time = 0;
    
    for (int l = 0; l <= max_level; l++)
    {
        for (auto obj = reverse_map [l].begin(); obj != reverse_map [l].end(); obj++)
        {
            naive_order [*obj] = time;
            time += 1;
        }
        
        if (l != max_level)
        {
            std::vector <Group*> next_level;
            
            for (auto obj = reverse_map [l].begin(); obj != reverse_map [l].end(); obj++)
            {
                vector <Group*> obj_kids;
                
                for (auto kid = reverse_map [l+1].begin(); kid != reverse_map [l+1].end(); kid++)
                {
                    if ((*obj)->nextGroups().find (*kid) != (*obj)->nextGroups().end ())
                        obj_kids.push_back (*kid);
                    
                    if (std::find (next_level.begin(), next_level.end(), *kid) == 
                        next_level.end ())
                    {
                        next_level.push_back (*kid);
                    }
                }
            }
            
            for (auto obj = reverse_map [l+1].begin(); obj != reverse_map[l+1].end(); obj++)
            {
                if (std::find (next_level.begin(), next_level.end(), *obj) == 
                    next_level.end ())
                {
                    next_level.push_back (*obj);
                }
            }
            
            reverse_map [l+1] = next_level;
        }
    }
    
    std::cout<< "optnaiveorder is"<<std::endl;
    for (auto it = naive_order.begin(); it != naive_order.end(); it++)
    {
        PyObject* func = PyObject_GetAttr (it->first->getCompAtIndex(0), str_func);
        PyObject* name = PyObject_GetAttr (func,
                                           Py_BuildValue ("s", "name"));
        PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
        const char* ss = PyBytes_AS_STRING (pyStr);
        std::cout<< "     "<< ss<< "  " << it->second << std::endl;
    }
}

string pystring_to_string (PyObject* str)
{
    char* _str;
    string cxx_str;
    
    PyObject* pyStr = PyUnicode_AsEncodedString(str, "utf-8", "Error ~");
    _str = PyBytes_AS_STRING (pyStr);
    //Py_DECREF (pyStr);
    cxx_str = string (_str);
    //free (_str);
    return cxx_str;
}

void classify_storage (const vector <PyObject*>& vec_comps, 
                       unordered_map <PyObject*, PyObject*>& comp_to_stg_class,
                       unordered_map <PyObject*, vector <PyObject*>>& new_storage_class_map)
{
    static PyObject* str_orig_stg_class = Py_BuildValue ("s", "orig_storage_class");
    
    //find_storage_equivalence
    unordered_map <string, vector <PyObject*>> storage_class_map;
    
    for (auto it = vec_comps.begin(); it != vec_comps.end(); it++)
    {
        string _strkey = stg_class_to_key [comp_to_orig_storage [*it]];
        if (storage_class_map.find (_strkey) == storage_class_map.end ())
        {
            storage_class_map [_strkey] = std::vector <PyObject*> ();
            storage_class_map [_strkey].push_back (*it);
        }
        else
        {
            storage_class_map [_strkey].push_back (*it);
        }
    }

    if (DEBUG)
    {
        std::cout<< "storage_class_map "<< std::endl;
        for (auto it = storage_class_map.begin(); it != storage_class_map.end (); it++)
        {
            std::cout<< "    "<<it->first << " " << it->second.size () << std::endl;
        }
    }
    
    //maximal_storage
    
    for (auto it = storage_class_map.begin(); it != storage_class_map.end (); it++)
    {
        static PyObject* str_func = Py_BuildValue ("s", "func");
        static PyObject* str_typ = Py_BuildValue ("s", "typ");
        static PyObject* str_ndims = Py_BuildValue ("s", "ndims");
        static PyObject* str_offsets = Py_BuildValue ("s", "offsets");
        static PyObject* str_generate_id = Py_BuildValue ("s", "generate_id");
        
        string key = it->first;
        vector <PyObject*>& class_comps = it->second;
        PyObject* helper_comp = class_comps [0];
        PyObject* func = comp_to_func [helper_comp];
        //check_and_print_exception ();
        PyObject* typ = func_to_typ [func];
        //check_and_print_exception ();
        long dims = func_to_ndims [func];
        //check_and_print_exception ();
        PyObject* helper_storage = comp_to_orig_storage[helper_comp];
        //check_and_print_exception ();
        vector <long> max_offset = comp_to_max_offsets [helper_comp];
        
        for (auto it2 = class_comps.begin(); it2 != class_comps.end(); it2++)
        {
            //PyObject* storage = comp_to_orig_storage[*it2];
            //check_and_print_exception ();
            //PyObject* offsets = PyObject_GetAttr (storage, str_offsets);
            //check_and_print_exception ();
            vector<long>& offset = comp_to_max_offsets [*it2];
            
            for (long i = 0; i < dims; i++)
            {
                long dim_off;
                
                dim_off = offset[i];
                //check_and_print_exception ();
                max_offset [i] = max (max_offset[i], dim_off);
            }
        }
        
        vector <PyObject*> vec_dim_sizes; //vector of tuples
        
        for (long dim = 0; dim < dims; dim++)
        {
            static PyObject* str_get_dim = Py_BuildValue ("s", "get_dim");
            static PyObject* str_coeff = Py_BuildValue ("s", "coeff");
            static PyObject* str_orig_param = Py_BuildValue ("s", "orig_param");
            PyObject* dim_storage;
            PyObject* get_dim;
            PyObject* new_size;
            PyObject* args;
            PyObject* max_offset_dim;
            
            max_offset_dim = PyLong_FromLong (max_offset [dim]);
            dim_storage = helper_storage_to_dim_storage [helper_storage][dim];
            args = PyTuple_Pack (2, dim_storage, max_offset_dim);
            //check_and_print_exception ();
            new_size = PyObject_CallObject (storage_mapping_get_dim_size, 
                                            args);
            Py_DECREF (args);
            Py_DECREF (max_offset_dim);
            ////std::cout<< "new_size " << new_size << std::endl;
            //check_and_print_exception ();
            
            vec_dim_sizes.push_back (PyTuple_Pack (2,
                                     PyObject_GetAttr (dim_storage, str_orig_param),
                                     new_size));
            //check_and_print_exception ();
        }
        
        PyObject* dim_sizes;
        
        dim_sizes = PyList_New (vec_dim_sizes.size ());
        //check_and_print_exception ();
    
        for (auto it2 = vec_dim_sizes.begin(); 
             it2 != vec_dim_sizes.end(); it2++)
        {
            PyList_SetItem (dim_sizes, (it2 - vec_dim_sizes.begin()), *it2);
            //check_and_print_exception ();
        }
        
        PyObject* max_storage;
        PyObject* py_dims;
        
        py_dims = PyLong_FromLong (dims);
        PyObject* args = PyTuple_Pack (3, typ, py_dims, dim_sizes);
        //check_and_print_exception ();
        max_storage = PyObject_CallObject (cls_Storage, args);
        //check_and_print_exception ();
        Py_DECREF (args);
        Py_DECREF (py_dims);
        //check_and_print_exception ();
        //args = PyTuple_Pack (0);
        //PyObject* generate_id;
        //generate_id = PyObject_GetAttr (max_storage, str_generate_id);
        //check_and_print_exception ();
        //PyObject_CallObject (generate_id, NULL);
        //check_and_print_exception ();
        new_storage_class_map [max_storage] = vector <PyObject*> ();
        
        for (auto it2 = class_comps.begin(); it2 != class_comps.end(); it2++)
        {
            comp_to_stg_class [*it2] = max_storage;
            new_storage_class_map [max_storage].push_back (*it2);
        }
        
        Py_DECREF (dim_sizes);
        
        for (auto it2 = vec_dim_sizes.begin(); 
             it2 != vec_dim_sizes.end (); it2++)
        {
            PyObject* tuple = *it2;
            PyObject* tuple1 = PyTuple_GetItem (tuple, 1);
            Py_DECREF (tuple);
            Py_DECREF (tuple1);
        }
    }
}

void getLivenessMap (const uint128_t hash_id, const vector<Group*>& vec_groups, 
                  unordered_map <Group*, int>& schedule, 
                  unordered_map <int, vector <PyObject*> >& liveness_map)
{
    for (auto it = schedule.begin(); it != schedule.end(); it++)
    {
        int last_live = -1;
        
        for (auto child = it->first->nextGroups().begin (); 
             child != it->first->nextGroups().end (); child++)
        {
            if (((*child)->hashID() & hash_id) == (*child)->hashID ())
            {
                last_live = max (last_live, schedule[*child]);
            }
        }
        
        if (liveness_map.find (last_live) == liveness_map.end ())
            liveness_map [last_live] = std::vector <PyObject*> ();
        
        liveness_map [last_live].push_back (it->first->getCompAtIndex(0));
    }
}

uint64_t getLiveoutsSize (uint128_t hash_id, vector <Group*> vec_groups, int& n_liveouts)
{
    uint64_t liveouts_size = 0;
    n_liveouts = 0;
    
    for (auto const& group: vec_groups)
    {
        for (auto const& next_group: group->nextGroups ())
        {
            if ((next_group->hashID () & hash_id) != next_group->hashID ())
            {
                PyObject* pygroup = next_group->getPyGroup ();
                uint64_t stg_size = 1;
                for (int i = 0; i < get_n_dimensions (pygroup); i++)
                {
                    //TODO: Optimize it because for every get_dim_reuse_for_pygroup call
                    //dict and list are accessed again. Make them access only once 
                    //and then iterate over list
                    stg_size *= get_dim_size_for_pygroup (pygroup, i);
                }
                
                liveouts_size += stg_size*IMAGE_ELEMENT_SIZE;
                n_liveouts++;
                break;
            }
        }
        
        if (group->nextGroups().size () == 0)
        {
            PyObject* pygroup = group->getPyGroup ();
            uint64_t stg_size = 1;
            
            for (int i = 0; i < get_n_dimensions (pygroup); i++)
            {
                //TODO: Optimize it because for every get_dim_reuse_for_pygroup call
                //dict and list are accessed again. Make them access only once 
                //and then iterate over list
                stg_size *= get_dim_size_for_pygroup (pygroup, i);
            }
                
            liveouts_size += stg_size*IMAGE_ELEMENT_SIZE;
            n_liveouts++;
        }
    }
    
    return liveouts_size;
}

uint64_t getTotalSizeUsed (uint128_t hash_id, int& number_of_buffers,
                           uint64_t& liveouts_size, uint64_t& livein_size)
{
    //TODO: Support TStencil 
    map <Group*, int, Group::GroupComparer> level_order;
    unordered_map <Group*, int> naive_order; //naive_order is schedule
    unordered_map <PyObject*, PyObject*> cmp_to_stg_class;
    unordered_map <PyObject*, vector <PyObject*>> storage_class_map;
    std::map <int, PyObject*> sorted_comps;
    unordered_map <PyObject*, vector<uint64_t> > array_pool; //Pool from computation to the vector of integers available
    uint128_t _hash_id = hash_id;
    vector <PyObject*> vec_comps;
    vector <Group*> vec_groups;
    vector <Group*> live_out_groups;
    unordered_set <Group*, Group::GroupHasher> live_in_groups;
    unordered_map <PyObject*, uint64_t> storage_map; //storage_map to store first index of array for each computation
    unordered_map <int, vector <PyObject*> > liveness_map;
    static PyObject* str_children = Py_BuildValue ("s", "children");
    static PyObject* str_func = Py_BuildValue ("s", "func");
    static PyObject* str_comps = Py_BuildValue ("s", "comps");
    uint64_t array_count = 0;
    int bit = 0;
    liveouts_size = 0;
    uint128_t nonLiveOutsHashID = 0;
    
    while (_hash_id != 0)
    {
        if (_hash_id & 1L == 1)
        {   
            uint128_t l = 1;
            Group* _g = Group::hashIDToGroup[l<<bit];
            PyObject* pyGroup = _g->getPyGroup ();
            
            if (!pyGroup)
            {
                _hash_id = _hash_id >> 1;
                bit++;
                continue;
            }
            
            if (!OptGroup::isLiveoutForGroup (hash_id, (OptGroup*)_g))
            {
                nonLiveOutsHashID |= (l << bit);
                vec_groups.push_back (_g);
                PyObject* comps = PyObject_GetAttr (pyGroup, str_comps);
                
                for (int j = 0; j < PyList_Size (comps); j++)
                {
                    PyObject* comp = PyList_GetItem (comps, j);

                    vec_comps.push_back (comp);
                }
            }
            else
                live_out_groups.push_back (_g);
            
            OptGroup::liveinsForChildInGroup (hash_id, (OptGroup*)_g, live_in_groups);
        }
        
        _hash_id = _hash_id >> 1;
        bit++;
    }
    
    get_level_order (nonLiveOutsHashID, vec_comps, vec_groups, level_order);
    naive_sched_objs (level_order, naive_order);
    if (DEBUG)
    {
        std::cout<<"optnaiveorder12121212: " << std::endl;
        for (auto it = naive_order.begin (); it != naive_order.end (); it++)
        {
            std::cout << it->second << " " << getPyCompFuncName (it->first->getCompAtIndex (0)) << std::endl;
        }
    }
    getLivenessMap (nonLiveOutsHashID, vec_groups, naive_order, liveness_map);
    if (DEBUG)
    {
        std::cout<<"optnaiveorderqqqqqqqqqq1212: " << std::endl;
        for (auto it = naive_order.begin (); it != naive_order.end (); it++)
        {
            std::cout << it->second << " " << getPyCompFuncName (it->first->getCompAtIndex (0)) << std::endl;
        }
    }
    classify_storage (vec_comps, cmp_to_stg_class, storage_class_map);
    liveouts_size = getLiveoutsSize (hash_id, live_out_groups, number_of_buffers);
    
    /*Get LiveIn Size*/
    for (auto const& __g: live_in_groups)
    {
        {
            PyObject* pygroup = __g->getPyGroup ();
            if (!pygroup)
                continue;
            uint64_t stg_size = 1;
            for (int i = 0; i < get_n_dimensions (pygroup); i++)
            {
                //TODO: Optimize it because for every get_dim_reuse_for_pygroup call
                //dict and list are accessed again. Make them access only once 
                //and then iterate over list
                stg_size *= get_dim_size_for_pygroup (pygroup, i);
            }
            
            livein_size += stg_size*IMAGE_ELEMENT_SIZE;
            break;
        }
    }
        
    if (DEBUG)
        std::cout << "liveout buffers "<< number_of_buffers << std::endl;
    /*//std::cout<< "optschedule is"<<std::endl;
    for (auto it = vec_comps.begin(); it != vec_comps.end(); it++)
    {
        PyObject* func = PyObject_GetAttr (*it, str_func);
        
        PyObject* name = PyObject_GetAttr (func,
                                           Py_BuildValue ("s", "name"));
        PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
        const char* ss = PyBytes_AS_STRING (pyStr);
            //std::cout<< "     "<< ss<< (it - vec_comps.begin()) << std::endl;
    }*/
    
    std::cout<< "cmp_to_stg_class is"<<std::endl;
    for (auto it = vec_comps.begin(); it != vec_comps.end(); it++)
    {
        PyObject* func = PyObject_GetAttr (*it, str_func);
        PyObject* name = PyObject_GetAttr (func,
                                           Py_BuildValue ("s", "name"));
        PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
        const char* ss = PyBytes_AS_STRING (pyStr);
        std::cout<< "     "<< ss<< (cmp_to_stg_class [*it]) << std::endl;
    }
    
    for (auto it = vec_comps.begin(); it != vec_comps.end(); it++)
    {
        array_pool [cmp_to_stg_class [*it]] = vector <uint64_t> ();
    }
    
    /*//std::cout<< "liveness map " << std::endl;
    
    for (auto it = liveness_map.begin(); it != liveness_map.end (); it++)
    {
        //std::cout<<it->first << std::endl;
        
        for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
        {
            PyObject* func = PyObject_GetAttr (*it2, str_func);
        
        PyObject* name = PyObject_GetAttr (func,
                                           Py_BuildValue ("s", "name"));
        PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
        const char* ss = PyBytes_AS_STRING (pyStr);
            //std::cout<< "     "<< ss<<std::endl;
        }
    }*/
    
    for (auto it = naive_order.begin(); it != naive_order.end(); it++)
    {
        sorted_comps[it->second] = it->first->getCompAtIndex(0);
    }
    
    if (DEBUG)
    {
        std::cout<<"sorted comps " << std::endl;
        for (auto it = sorted_comps.begin (); it != sorted_comps.end (); it++)
        {
            std::cout << it->first << " " << getPyCompFuncName (it->second) << std::endl;
        }
        
        std::cout<<"optnaiveorder: " << std::endl;
        for (auto it = naive_order.begin (); it != naive_order.end (); it++)
        {
            std::cout << it->second << " " << getPyCompFuncName (it->first->getCompAtIndex (0)) << std::endl;
        }
    }
    
    unordered_map <PyObject*, int> n_stg_class_arrays; //Number of arrays for each storage class
    
    for (auto it = sorted_comps.begin(); it != sorted_comps.end(); it++)
    {
        PyObject* comp = it->second;
        PyObject* storage_class = cmp_to_stg_class [comp];
        PyObject* func = comp_to_func [comp];
        
        uint64_t num_reqd, num_available, num_allocated, deficit;
        bool is_tstencil_grp = false;
        vector <uint64_t> allocated_arrays;
        int time;
        
        if (PyObject_IsInstance (func, tstencil_cls))
        {
            is_tstencil_grp = true;
        }
        
        if (is_tstencil_grp == true)
            num_reqd = 2;
        else
            num_reqd = 1;
        
        if (DEBUG)
        {
            std::cout << "comp is " << getPyCompFuncName (comp) << " sched is " << it->first << std::endl;
            std::cout<< "num_reqd " << num_reqd << std::endl;
        }
        num_available = array_pool [storage_class].size ();
        if (DEBUG)
            std::cout<< "num available " << num_available <<std::endl;
        num_allocated = min (num_available, num_reqd);
        
        for (uint64_t i = 0; i < num_allocated; i++)
        {
            allocated_arrays.push_back (array_pool[storage_class].back ());
            array_pool[storage_class].pop_back ();
        }
        
        deficit = num_reqd - num_allocated;
        
        if (deficit > 0)
        {
            for (uint64_t i = 0; i < deficit; i++)
            {
                allocated_arrays.push_back (array_count + 1 + i);
            }
            
            array_count += deficit;
            n_stg_class_arrays [storage_class] += deficit;
        }
        
        storage_map [comp] = allocated_arrays[0];
        time = it->first;
        
        if (liveness_map.find (time) != liveness_map.end ())
        {
            ////std::cout<< "time " << time << " in liveness_map " << std::endl;
            vector <PyObject*> free_comps;
            free_comps = liveness_map [time];
            PyObject* cmp_stg_class;
            
            for (auto it2 = free_comps.begin(); it2 != free_comps.end(); it2++)
            {
                cmp_stg_class = cmp_to_stg_class[*it2];
                /*PyObject* func = PyObject_GetAttr (*it2, str_func);
        
        PyObject* name = PyObject_GetAttr (func,
                                           Py_BuildValue ("s", "name"));
        PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
        const char* ss = PyBytes_AS_STRING (pyStr);
            //std::cout<< "free comp "<< ss<<std::endl;
            */
                array_pool [cmp_stg_class].push_back (storage_map [*it2]);
            }
        }
    }

    std::cout<< "optarray pool" << std::endl;
    for (auto it = array_pool.begin(); it != array_pool.end(); it++)
    {
        PyObject* name = PyObject_Str (it->first);
        PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
        const char* ss = PyBytes_AS_STRING (pyStr);
        std::cout<<ss<<std::endl;
        
        for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
        {
            std::cout<< "   " << *it2 << std::endl;
        }
    }
    
    if (DEBUG)
            std::cout<< "Array count " << array_count << " " << std::bitset<41>((uint64_t) hash_id) << std::endl;
    uint64_t total_size = 0;
    if (DEBUG)
            std::cout<< "total size " << n_stg_class_arrays.size () << std::endl;
    for (auto it = n_stg_class_arrays.begin(); 
         it != n_stg_class_arrays.end (); it++)
    {
        PyObject* stg_class = it->first;
        int n_arrays = it->second;
        if (DEBUG)
            std::cout << "n_arrays " << n_arrays << std::endl;
        PyObject* pygroup = pyGroupForComp [storage_class_map [stg_class][0]];
        if (DEBUG)
            std::cout << "pygroup " << pygroup << std::endl;
        uint64_t stg_size = 1;
        std::cout << "dim size : ";
        for (int i = 0; i < get_n_dimensions (pygroup); i++)
        {
            //TODO: Optimize it because for every get_dim_reuse_for_pygroup call
            //dict and list are accessed again. Make them access only once 
            //and then iterate over list
            uint64_t __size = get_dim_size_for_pygroup (pygroup, i);
            std::cout << "(dim: " << i << ", size: " << __size << ")" << std::endl;
            stg_size *= __size;
        }
        std::cout << std::endl;
        total_size += stg_size*n_arrays*IMAGE_ELEMENT_SIZE;
    }
    
    unordered_set <PyObject*> stg_class_removed;
    
    for (auto it = storage_class_map.begin(); it != storage_class_map.end(); it++)
    {
        Py_DECREF (it->first);
    }
    
    number_of_buffers += array_count;
    
    return total_size + liveouts_size;
}

inline uint64_t dim_size_std_dev (std::vector <std::vector <uint64_t> >& dim_size_diff, int max_dim)
{
    std::vector <double> sum, mean; //vectors for sum and mean in each dimension
    sum = std::vector <double> (max_dim, 0);
    mean = std::vector <double> (max_dim, 0);
    
    for (auto it_comp = dim_size_diff.begin (); it_comp != dim_size_diff.end (); 
         it_comp++)
    {
        for (int dim = 0; dim < (*it_comp).size (); dim++)
        {
            sum[dim] += (*it_comp)[dim];
        }
    }
    
    for (int i = 0; i < sum.size (); i++)
    {
        mean[i] = sum[i]/dim_size_diff.size (); //sum [i]->sum of dim sizes in each dimension
        //dim_size_diff.size () -> number of computations
    }
    
    //Calculate mean difference in each dimension
    std::vector <double> mean_diff = std::vector <double> (max_dim, 0);
    
    for (auto it_comp = dim_size_diff.begin (); it_comp != dim_size_diff.end (); 
         it_comp++)
    {
        for (int dim = 0; dim < (*it_comp).size (); dim++)
        {
            mean_diff [dim] += abs ((*it_comp)[dim]- mean[dim])/dim_size_diff.size ();///mean[dim]*500.0f;
        }
    }
    
    //return sum of mean diff
    //return ceil (*(mean_diff.end()-1)*100.0f/(*(mean.end()-1)) + 
    //    ((std::accumulate (mean_diff.begin(), mean_diff.end ()-1, 0.0)*100.0f) / std::accumulate (mean.begin(), mean.end ()-1, 0.0)));
    //std::accumulate (mean_diff.begin(), mean_diff.end (), 0.0); for campipe
    double sum_mean_dim = std::accumulate (mean.begin(), mean.end (), 0.0);
    double sum_mean_diff = std::accumulate (mean_diff.begin(), mean_diff.end (), 0.0);
    std::cout <<"mean sum" << std::accumulate (mean.begin(), mean.end (), 0.0) << " mean_diff sum " <<std::accumulate (mean_diff.begin(), mean_diff.end (), 0.0) <<std::endl;
    if (sum_mean_diff/sum_mean_dim > 0.1f)
        return sum_mean_diff;//\ /std::accumulate (mean.begin(), mean.end (), 0.0);
    else
      return (sum_mean_diff*100.0f)/sum_mean_dim;
}

int getNumberOfEdges (uint128_t _hash_id)
{
    vector <Group*> vec_groups;
    int n_edges  = 0;
    int bit = 0;
    
    while (_hash_id != 0)
    {
        if (_hash_id & 1L == 1)
        {   
            uint128_t l = 1;
            Group* _g = Group::hashIDToGroup[l<<bit];
            PyObject* pyGroup = _g->getPyGroup ();
            
            if (!pyGroup)
            {
                _hash_id = _hash_id >> 1;
                bit++;
                continue;
            }
            
            for (auto const &next_group : _g->nextGroups ())
            {
                if ((next_group->hashID () & _hash_id) == next_group->hashID ())
                    n_edges++;
            }
        }
        
        _hash_id = _hash_id >> 1;
        bit++;
    }
    
    return n_edges;
}

inline uint64_t cost (uint128_t hash_id, std::vector <uint64_t>& tile_sizes)
{
    uint64_t bit = 0;
    uint64_t nOnes = 0;
    int total_comps = 0;
    bool is_reduction = false;
    bool is_small_grp = false;
    bool is_tstencil_grp = false;
    bool is_const_grp = false;
    bool is_dummy_source_group = false;
    uint128_t _hash_id = hash_id;
    static PyObject* str_comps = Py_BuildValue ("s", "comps");
    static PyObject* str_func = Py_BuildValue ("s", "func");
    static PyObject* str_children = Py_BuildValue ("s", "children");
    static PyObject* str_is_const_func = Py_BuildValue ("s", "is_const_func");
    PyObject* list_groups_for_overlap = PyList_New (0);
    PyObject* list_overlap_size;
    int max_dim = get_max_dimensions_for_group (hash_id);
    int max_dims = max_dim;
    uint64_t max_dim_reuse = 0;
    uint64_t tile_sizes_product = 0;
    uint64_t max_dim_tile_size;
    int n_zero_dim_reuse = 0;
    int n_non_zero_dims = 0;
    std::vector <uint64_t> dim_reuse = std::vector <uint64_t> (max_dims, 0);
    std::vector <uint64_t> dim_size = std::vector <uint64_t> (max_dims, 0);
    int n_threads = 0, liveouts_buffer= 0;
    int _cost;
    int n_buffers = 0;
    bool all_n_dims_equal = true;
    int last_group_dim = -1;
    uint64_t liveouts_size = 0;
    uint64_t liveins_size = 0;
    std::vector <std::vector <uint64_t> > dim_size_diff; // 2-D vector with rows as number of comps and cols as max_dim
    uint64_t tile_size = 0;
    uint64_t totalsizeused = getTotalSizeUsed (hash_id, n_buffers, 
                                               liveouts_size, liveins_size);
    
    if (DEBUG)
        std::cout << "n_buffers " << n_buffers << std::endl;
    if (DEBUG)
        std::cout<< "totalsizeused " << totalsizeused << std::endl;
    
    if (totalsizeused == 0 && n_buffers == 0)
    {
        //Dummy group
        return 0;
    }
    
    while (_hash_id != 0)
    {
        if (_hash_id & 1L == 1)
        {
            uint128_t l = 1;
            PyObject* pyGroup = Group::hashIDToGroup[l<<bit]->getPyGroup ();
            nOnes++;
            if (pyGroup != NULL)
            {                
                
                PyObject* comps = PyObject_GetAttr (pyGroup, str_comps);
                
                total_comps += PyList_Size (comps);
                
                for (int j = 0; j < PyList_Size (comps); j++)
                {
                    PyObject* comp = PyList_GetItem (comps, j);
                    PyObject* func = PyObject_GetAttr (comp, str_func);
                    if (PyObject_IsInstance (func,
                                             reduction_cls))
                    {
                        is_reduction = true;
                        
                    }
                    
                    if (PyObject_IsInstance (func,
                                             tstencil_cls))
                    {
                        is_tstencil_grp = true;
                    }
                       
                    PyObject* is_const_func = PyObject_GetAttr (func, str_is_const_func);
                    if (is_const_func == Py_True)
                    {
                        is_const_grp = true;
                    }
                       
                    if (small_comps_set.find (comp) != small_comps_set.end())
                    {
                        is_small_grp = true;
                    }
                }
                
                dim_size_diff.push_back (std::vector <uint64_t> (max_dim, 0));
                
                for (int i = 0; i < get_n_dimensions (pyGroup); i++)
                {
                    //TODO: Optimize it because for every get_dim_reuse_for_pygroup call
                    //dict and list are accessed again. Make them access only once 
                    //and then iterate over list
                    //PyObject* comp = PyList_GetItem (comps, 0);
                    //PyObject* func = PyObject_GetAttr (comp, str_func);
                    dim_reuse [i] += get_dim_reuse_for_pygroup (pyGroup, i)*IMAGE_ELEMENT_SIZE;
                    uint64_t _dim_size = get_dim_size_for_pygroup (pyGroup, i)*IMAGE_ELEMENT_SIZE;
                    dim_size [i] += _dim_size;
                    dim_size_diff [dim_size_diff.size () -1][i] = _dim_size;
                }
                
                if (all_n_dims_equal == true)
                {
                    if (last_group_dim != -1 and last_group_dim != get_n_dimensions (pyGroup))
                    {
                        all_n_dims_equal = false;
                    }
                    
                    last_group_dim = get_n_dimensions (pyGroup);
                }
                
                PyList_Append (list_groups_for_overlap, pyGroup);
            }
            else
            {
                is_dummy_source_group = true;
            }
        }
        
        _hash_id = _hash_id >> 1;
        bit++;
    }
    
    _hash_id = hash_id;
    
    if ((is_reduction || is_const_grp || is_small_grp || is_tstencil_grp || is_dummy_source_group) && nOnes == 1)
        return 1L;
    
    if (is_reduction || is_const_grp || is_small_grp || is_tstencil_grp || is_dummy_source_group)
        return -1;
    
    if (total_comps > grp_size)
        return -1;
    
    if (is_dummy_source_group)
        return -1;
        
    if (nOnes == 1)
        return INT_MAX;
    
    
    #define INCLUDE_OVERLAP
    #if INCLUDE_OVERLAP
    PyObject* pytotalsizeused = PyLong_FromLong (totalsizeused);
    PyObject* pyn_buffers = PyLong_FromLong (n_buffers);
    PyObject* args = PyTuple_Pack (3, list_groups_for_overlap, pytotalsizeused, pyn_buffers);
    PyObject* tuple_return = PyObject_CallObject (get_overlapping_size_func, args);
    check_and_print_exception ();
    PyObject* overlap_obj = PyTuple_GetItem (tuple_return, 0);
    check_and_print_exception ();
    PyObject* pytile_size = PyTuple_GetItem (tuple_return, 1);
    check_and_print_exception ();
    uint64_t overlap_size = PyLong_AsLong (overlap_obj)*IMAGE_ELEMENT_SIZE;
    tile_size = PyLong_AsLong(pytile_size)*IMAGE_ELEMENT_SIZE;
    if (tile_size == 0 && overlap_size == -4)
    {
        return -1;
    }
    
    if (tile_size == 0 && overlap_size == 0)
    {
        tile_size = totalsizeused/N_CORES;
        if (tile_size > L2_CACHE_SIZE)
            tile_size = L2_CACHE_SIZE;
        
        overlap_size = 0;
    }
    
    if (DEBUG)
        std::cout << "tile_size "<< tile_size << std::endl;
    assert (tile_size != 0);
    Py_DECREF (args);
    //Py_DECREF (tuple_return);
    Py_DECREF (pytotalsizeused);
    Py_DECREF (list_groups_for_overlap);
    Py_DECREF (pyn_buffers);
    #else
    uint64_t overlap_size = 0;
    #endif
    
     if (DEBUG)
            std::cout<< "totalsizeused " << totalsizeused << " tile_size for " << std::bitset<41>((uint64_t)_hash_id)<<" " << tile_size << " liveins_size " << liveins_size <<std::endl;
    
    uint64_t mean_dim_diff = dim_size_std_dev (dim_size_diff, max_dims);

    if (L2_CACHE_SIZE < totalsizeused/N_CORES)
    {
        n_threads = totalsizeused/L2_CACHE_SIZE;
        
        if (DEBUG)
            std::cout<< "L2 threads for " << std::bitset<41>((uint64_t)_hash_id)<<" " << n_threads << std::endl;
    }
    else
        n_threads = N_CORES;

    if (DEBUG)    
        std::cout << "_cost mean_dim_diff " << mean_dim_diff << " live ratio " << ((liveins_size+liveouts_size)*100)/(tile_size*n_buffers+liveins_size) <<
            " threads ratio " << (((n_threads+N_CORES-1)%N_CORES)*100)/N_CORES << " all_n_dims_equal " << all_n_dims_equal << " overlap " << overlap_size << " tile_size " << tile_size << std::endl;

    _cost = mean_dim_diff + ((liveins_size+liveouts_size))/(tile_size*n_buffers) - 
            (((n_threads+N_CORES-1)%N_CORES)*100.0f)/N_CORES;
            
    if (((100*overlap_size)/tile_size) <= 150)
    {
        /*
         * TODO:
         * If above condition is not set then the cost due to
         * overlap size of main grouping of Harris becomes very high
         * and that grouping is not returned. Should be easy to work 
         * around by changing weights.
         * */
        _cost += (1000.0f*50.0f*1.5*overlap_size)/(tile_size*n_buffers);
    }
    
    if (!all_n_dims_equal)
        _cost += 500;
        
    if (DEBUG)
        std::cout << "total cost " << _cost << std::endl;
    
    if (overlap_size > tile_size)
        ;//_cost += 100000;
    
    if (total_comps <= grp_size && (overlap_size < 1.5*(tile_size) 
        || (overlap_size > 1.7*(tile_size) && overlap_size < 2*(tile_size))))
    {
        /*
         * TODO:
         * A Little Hack here, < 1.5*(tile_size) helps with pyramid_blend and harris
         * And interval of (1.7, 2)*tile_size helps with bilateral_grid.
         * Basically we do not have a single (1, n)*tile_size range where
         * these benchmarks will run with the best grouping.
         * Should be easy to work around by changing weights.
         */
       return _cost ;
    }
    
    //Invalid Grouping
    return -1;
}

inline uint64_t cost2 (uint128_t hash_id, std::vector <uint64_t>& tile_sizes)
{
    uint64_t bit = 0;
    uint64_t nOnes = 0;
    long totalComps = 0;
    bool is_reduction = false;
    bool is_small_grp = true;
    bool is_tstencil_grp = false;
    bool is_const_grp = false;
    uint128_t _hash_id = hash_id;
    static PyObject* str_comps = Py_BuildValue ("s", "comps");
    static PyObject* str_func = Py_BuildValue ("s", "func");
    static PyObject* str_children = Py_BuildValue ("s", "children");
    static PyObject* str_is_const_func = Py_BuildValue ("s", "is_const_func");

    while (hash_id != 0)
    {
        if (hash_id & 1L == 1)
        {   
            nOnes++;
            uint128_t l = 1;
            PyObject* pyGroup = Group::hashIDToGroup[l<<bit]->getPyGroup ();
            //std::cout<< "pygroup " << pyGroup << std::endl;
            if (pyGroup != NULL)
            {    
                PyObject* list = PyObject_GetAttr (pyGroup, str_comps);
                PyObject* name = PyObject_GetAttr (pyGroup,
                                                   Py_BuildValue ("s", 
                                                                  "name"));
                PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
                const char* ss = PyBytes_AS_STRING (pyStr);
                //std::cout<<ss << " ";
                totalComps += PyList_Size (list);
                
                PyObject* comps = PyObject_GetAttr (pyGroup, str_comps);
                    //if (!comps)
                    //    continue;
                        
                    for (int j = 0; j < PyList_Size (comps); j++)
                    {
                       PyObject* comp = PyList_GetItem (comps, j);
                       PyObject* func = PyObject_GetAttr (comp, str_func);
                       if (PyObject_IsInstance (func,
                                                reduction_cls))
                       {
                            is_reduction = true;
                            //PyObject* name = PyObject_GetAttr (func,
                            //                       Py_BuildValue ("s", 
                            //                                      "name"));
                //PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
                //const char* ss = PyBytes_AS_STRING (pyStr);
                //////std::cout<<"REDUCTIONwe "<<ss << " ";
                       }
                       if (PyObject_IsInstance (func,
                                                tstencil_cls))
                       {
                            is_tstencil_grp = true;
                       }
                       
                       PyObject* is_const_func = PyObject_GetAttr (func, str_is_const_func);
                       if (is_const_func == Py_True)
                       {
                            is_const_grp = true;
                       }
                       
                       if (small_comps_set.find (comp) == small_comps_set.end())
                       {
                            is_small_grp = false;
                       }
                    }
                
                PyObject* children = PyObject_GetAttr (pyGroup, str_children);
                //if (!children)
                //    continue;
                    
                for (int i = 0; i < PyList_Size (children); i++)
                {
                    PyObject* child = PyList_GetItem (children, i);
                    //if (!child)
                    //    continue;
                        
                    PyObject* comps = PyObject_GetAttr (child, str_comps);
                    //if (!comps)
                    //    continue;
                        
                    for (int j = 0; j < PyList_Size (comps); j++)
                    {
                       PyObject* comp = PyList_GetItem (comps, j);
                       PyObject* func = PyObject_GetAttr (comp, str_func);
                       if (PyObject_IsInstance (func,
                                                reduction_cls))
                       {
                            is_reduction = true;
                       }
                       if (PyObject_IsInstance (func,
                                                tstencil_cls))
                       {
                            is_tstencil_grp = true;
                       }
                       
                       PyObject* is_const_func = PyObject_GetAttr (func, str_is_const_func);
                       if (is_const_func == Py_True)
                       {
                            is_const_grp = true;
                       }                       
                    }
                }
            }
        }
        
        hash_id = hash_id >> 1;
        bit++;
    }
    //std::cout<< "cost is here " << std::endl;
    
    if (totalComps <= grp_size && !is_reduction && !is_const_grp && !is_small_grp && !is_tstencil_grp){////std::cout<<"hash_id9 " << std::bitset<8>((int)_hash_id) << " reduction " << is_reduction<< " nOnes " << nOnes << std::endl;
        return 1L;}
    //else
        //////std::cout<< "totalComps " << totalComps << std::endl;
    ////std::cout<<"hash_id " << std::bitset<8>((int)_hash_id) << " reduction " << is_reduction<< " nOnes " << nOnes << std::endl;
    if (is_reduction && nOnes == 1) //Group contains only one group, which happens to be the reduction too
        return 1L;
    ////std::cout<<"hash_id23 " << std::bitset<8>((int)_hash_id) << " reduction " << is_reduction<< " nOnes " << nOnes << std::endl;
    return -1;
}

inline bool reachable (uint128_t start, uint128_t end, bool isEndGroup)
{
    queue<uint128_t> q;
    ////std::cout<< "reachable isEndGroup " << isEndGroup << " start " << std::bitset<41>((uint64_t)start) << " end " << std::bitset<41>((uint64_t)end)<<std::endl;
    q.push (start);
    
    while (!q.empty ())
    {
        //std::cout<< "q.size () "<< q.size () << std::endl;
        vector<uint128_t>* n;
        uint128_t hash_id = q.front ();
        q.pop();
        //std::cout<< "popped " << std::bitset<41>((uint64_t)hash_id) << std::endl;
        
        if (isEndGroup)
        {
            if ((hash_id & end) == hash_id)
               return true;
        }
        else
            if (hash_id == end)
                return true;
        
        n = OptGroup::nextGroupsHashID (hash_id);
        
        for (auto it = n->begin(); it != n->end(); it++)
        {
            q.push (*it);
        }
        
        delete n;
    }

    return false;
}

inline bool isCycle (uint128_t hash_id, uint128_t g, vector<uint128_t>* next)
{
    //std::cout<< "isCycle next NOT null hash_id "<< std::bitset<41>((uint64_t)hash_id) << " g " << std::bitset<41>((uint64_t)g)<<std::endl;
    
    for (auto it = next->begin (); it != next->end(); it++)
    {
        if (*it != g && reachable (*it, g, false))
            return true;
    }
    
    return false;
}

inline bool isCycle (uint128_t hash_id, uint128_t g)
{
    vector<uint128_t>* next = OptGroup::nextGroupsHashID(g);
    //std::cout<< "isCycle next null hash_id "<< std::bitset<41>((uint64_t)hash_id) << " g " << std::bitset<41>((uint64_t)g)<<std::endl;
    for (auto it = next->begin (); it != next->end(); it++)
    {
        if (reachable (*it, hash_id, true))
            return true;
    }
    
    return false;
}

void getAllNodesInNextPath (uint128_t node, uint128_t dest, 
                            unordered_set<uint128_t, uint128Hasher>& nodes_in_path,
                            vector<uint128_t>& path,
                            vector<uint128_t>* (*next_groups) (uint128_t))
{
    path.push_back (node);
    
    if (node == dest)
    {
        for (auto it = path.begin (); it != path.end (); it++)
        {
            nodes_in_path.insert (*it);
        }
    }
    else
    {
        vector<uint128_t>* n;
        
        n = next_groups (node);
        
        for (auto it = n->begin(); it != n->end(); it++)
        {
            if (nodes_in_path.find (*it) == 
                nodes_in_path.end())
                getAllNodesInNextPath (*it, dest, nodes_in_path, path, next_groups);
        }
        
        delete n;
    }
    
    path.pop_back ();
}

inline uint128_t getCycleGroup (uint128_t curr_group, uint128_t next_group)
{
    stack<uint128_t> _stack;
    vector<uint128_t> path_vector;
    unordered_set<uint128_t, uint128Hasher> visited_curr_to_next;
    unordered_set<uint128_t, uint128Hasher> visited_next_to_curr;
    vector<uint128_t> intersection;
    uint128_t cg = curr_group | next_group;
    vector<uint128_t>* next_hashids = OptGroup::nextGroupsHashID (curr_group);
    vector<uint128_t>* prev_hashids;
    
    for (auto it = next_hashids->begin(); it != next_hashids->end(); it++)
    {
        if (*it != next_group)
            //_stack.push (*it);
            getAllNodesInNextPath (*it, next_group, visited_curr_to_next, 
                                   path_vector, OptGroup::nextGroupsHashID);
    }
    
    delete next_hashids;

    /*while (!_stack.empty ())
    {
        //std::cout<< "_stack.size () "<< _stack.size () << std::endl;
        vector<uint128_t>* n;
        uint128_t hash_id = _stack.top ();
        _stack.pop();
        //std::cout<< "popped " << std::bitset<41>((uint64_t)hash_id) << std::endl;
        
        if (hash_id == next_group)
        {
            for (auto it = path_vector.begin(); it != path_vector.end();
                 it++)
                visited_curr_to_next.insert (*it);
                
            path_vector.pop_back();
        }
        
        path_vector.push_back (hash_id);
        
        n = OptGroup::nextGroupsHashID (hash_id);
        
        for (auto it = n->begin(); it != n->end(); it++)
        {
            if (visited_curr_to_next.find (*it) == 
                visited_curr_to_next.end())
                _stack.push (*it);
        }
        
        delete n;
    }*/
    
    
    
    /*while (!_stack.empty ())
        _stack.pop();
    
    path_vector.clear ();
    prev_hashids = OptGroup::prevGroupsHashID (next_group);
    
    for (auto it = prev_hashids->begin(); it != prev_hashids->end(); it++)
    {
        if (*it != curr_group)
            getAllNodesInNextPath (*it, curr_group, visited_next_to_curr, 
                                   path_vector, OptGroup::prevGroupsHashID);
    }
    
    delete prev_hashids;

    while (!_stack.empty ())
    {
        //std::cout<< "prev _stack.size () "<< _stack.size () << std::endl;
        vector<uint128_t>* n;
        uint128_t hash_id = _stack.top ();
        visited_next_to_curr.insert (hash_id);
        _stack.pop();
        //std::cout<< "prevpopped " << std::bitset<41>((uint64_t)hash_id) << std::endl;
        
        if (hash_id == curr_group)
        {
            //std::cout<< "hash_id == curr_group " << std::endl;
            for (auto it = path_vector.begin(); it != path_vector.end();
                 it++)
            {
                //std::cout<< "prevpopped " << std::bitset<41>((uint64_t)hash_id) << std::endl;
                visited_next_to_curr.insert (*it);
            }
            //std::cout<< "added to visited next to curr " << std::endl;    
            path_vector.pop_back();
        }
        
        path_vector.push_back (hash_id);
        
        n = OptGroup::prevGroupsHashID (hash_id);
        //std::cout<< "343434343sr " << std::endl;
        for (auto it = n->begin(); it != n->end(); it++)
        {
            if (visited_next_to_curr.find (*it) == 
                visited_next_to_curr.end())
                _stack.push (*it);
        }
        
        delete n;
    }*/
    
    //std::cout<< "visited_curr_to_next " << std::endl;
    
    for (auto it = visited_curr_to_next.begin (); it != visited_curr_to_next.end (); it++)
    {
        //std::cout<< "   " << std::bitset<41>((uint64_t)*it) << std::endl;
    }
    
    /*//std::cout<< "visited_next_to_curr " << std::endl;
    
    for (auto it = visited_next_to_curr.begin (); it != visited_next_to_curr.end (); it++)
    {
        //std::cout<< "   " << std::bitset<41>((uint64_t)*it) << std::endl;
    }
    
    for (auto it = visited_curr_to_next.begin (); it != visited_curr_to_next.end (); it++)
    {
        if (visited_next_to_curr.find (*it) != visited_next_to_curr.end())
            intersection.push_back (*it);
    }
    
    //std::cout<< "intersection " << std::endl;
    
    for (auto it = intersection.begin (); it != intersection.end (); it++)
    {
        //std::cout<< "   " << std::bitset<41>((uint64_t)*it) << std::endl;
    }*/
    
    for (auto it = visited_curr_to_next.begin(); it != visited_curr_to_next.end(); it++)
    {
        cg |= *it;
    }
    
    return cg;
}

inline uint64_t cfgMemoize (uint128_t hash_id)
{
    if (DEBUG)
        std::cout << "STARTING for hash_id " << std::bitset<41> ((uint64_t)hash_id) << std::endl;
    vector<uint128_t>* n;
    
    auto g = T.find (hash_id);
    if (g != T.end())
    {
        ////std::cout<<"found hash_id " << std::bitset<41>((uint64_t)hash_id)<<std::endl;
        foundInDict += 1;
        return g->second;
    }
    
    //std::cout<<"to get next groups "<< std::bitset<41>((uint64_t)hash_id)<<std::endl;
    n = OptGroup::nextGroupsHashID (hash_id);
    
    /*////std::cout<< "next_groups: ";
    for (auto it = n->begin(); it != n->end (); ++it)
        ////std::cout<< (*it) << " ";
    ////std::cout<<std::endl;*/
    //std::cout<<"hash_id " << std::bitset<41>((uint64_t)hash_id)<<std::endl;
    if (n->size() == 0)
    {
        //T[hash_id] = cost(hash_id);
        ////std::cout<<" co2222st " <<cost(hash_id) << std::endl;
        delete n;
        int max_dim = get_max_dimensions_for_group (hash_id);
        std::vector <uint64_t> tile_sizes = std::vector <uint64_t> (max_dim, 0);
        hashIDToOptHashID[hash_id] = hash_id;
        uint64_t _cost = cost(hash_id, tile_sizes);
        optHashIDToTileSizes [hash_id] = tile_sizes;
        
        /*//std::cout<< "SETTING121 TILE SIZES FOR " << std::bitset<41>((uint64_t)hash_id)<< std::endl;
        for (auto it = tile_sizes.begin (); it != tile_sizes.end (); it++)
        {
            //std::cout<< "   " << *it << std::endl;
        }*/
        
        if (hash_id && (!(hash_id & (hash_id-1))))
            return _cost;
            
        return _cost;
        //return 0;
    }
    
    uint64_t totalMinCost;
    
    runningTime += 1;
    
    //if ((runningTime & 16283) == 0)
    {
        ////std::cout<< runningTime << " " << foundInDict<<std::endl;
    }
    
    if (numberOfOnes (hash_id) < (1 << logMaxChildren))
    {
        totalMinCost = 1L << 60;
        uint128_t optHashID = 0;
        uint128_t next_hash_id = 0;
        std::vector <uint64_t> opt_tile_sizes;
        std::vector <uint64_t> opt_tile_sizes_next_group;
        std::vector <uint64_t> tile_sizes_curr_group;
        int max_dim_curr_group;
        uint64_t cost_curr_group;
        uint64_t cost_next_groups = 0;
        unordered_set <uint128_t, uint128Hasher> includedVertices;
        
        for (auto it = n->begin(); it != n->end(); ++it)
        {
            if (includedVertices.find (*it) == includedVertices.end())
            {
                std::vector <uint64_t> tile_sizes_next_group;
                int max_dim_next_group = get_max_dimensions_for_group (*it);
                tile_sizes_next_group = std::vector <uint64_t> (max_dim_next_group, 0);
                //std::cout<<"filled it " << std::bitset<41>((uint64_t)*it)<<std::endl;
                uint64_t q = cfgMemoize (*it);
                cost_next_groups += q;
                ////std::cout<<"filled it " << std::bitset<41>((uint64_t)*it)<<" cost "<< q + cost(*it)<<std::endl;     
                if (DEBUG)
                    std::cout << "next " << std::bitset<41> ((uint64_t)(*it)) << " cost " << q << std::endl;
                
                //includedVertices.insert (*it);
            }
        }
        
        
        if (DEBUG)
            std::cout << "cost_next_groups " << cost_next_groups << std::endl;
        max_dim_curr_group = get_max_dimensions_for_group (hash_id);
        tile_sizes_curr_group = std::vector <uint64_t> (max_dim_curr_group, 0);
        cost_curr_group = cost(hash_id, tile_sizes_curr_group);
        //std::cout<< "serw33434 " << cost_curr_group << std::endl;
        //if (hash_id == 0b00000000000000000000000000000000011000000)
            ;//cost_next_groups -= 833609;
        for (auto it = n->begin(); it != n->end(); ++it)
        {
            std::vector <uint64_t> tile_sizes_cg;
            std::vector <uint64_t> tile_sizes_next_group;
            
            int max_dim_next_group;
            
            uint128_t cg = (*it) | hash_id;
            if (DEBUG)
                std::cout<<"it " << std::bitset<41>((uint64_t)*it)<<std::endl;
            uint64_t cost1;
            if (isCycle (hash_id, *it, n))
            {
                //std::cout<< "cycle found in grouping "<<std::endl;
                /*cg = getCycleGroup (hash_id, *it);
                //std::cout<< "new cg is " << std::bitset <41>((uint64_t)cg)<<std::endl;
                
                int max_dim_cg = get_max_dimensions_for_group (cg);
                //std::cout<< "max_dim_cg " << max_dim_cg << std::endl;
                tile_sizes_cg = std::vector <uint64_t> (max_dim_cg, 0);
                
                if (cost (cg, tile_sizes_cg) == -1)
                    cost1 = -1;
                else
                {
                    //std::cout<< "memoizing " << std::endl;
                    uint64_t o = cfgMemoize (cg);
                    cost1 = o;
                }*/
                cost1 = -1;
            }
            else if (isCycle (hash_id, *it))
            {
                //std::cout<< "cycle found in grouping2 "<<std::endl;
                cg = getCycleGroup (hash_id, *it);
                //std::cout<< "new cg is2 " << std::bitset <41>((uint64_t)cg)<<std::endl;
                cost1 = -1; //Invalid Grouping
            }
            else
            {   
                //std::cout<< "getting cost for cg " << std::endl;
                int max_dim_cg = get_max_dimensions_for_group (cg);
                //std::cout<< "max_dim_cg " << max_dim_cg << std::endl;
                tile_sizes_cg = std::vector <uint64_t> (max_dim_cg, 0);
                uint64_t cg_cost = cost (cg, tile_sizes_cg);
                //std::cout<< "cg_cost " << (int64_t)cg_cost << std::endl;
                if ((int64_t)cg_cost == -1)
                    cost1 = -1;
                else
                {
                    if (DEBUG)
                        std::cout<< "memoizing " << std::endl;
                    uint64_t o = cfgMemoize (cg);
                    cost1 = o;
                }
            }
            
            //std::cout<<"sdfsdfsd44443535 "<< std::endl;
            max_dim_next_group = get_max_dimensions_for_group (*it);
            tile_sizes_next_group = std::vector <uint64_t> (max_dim_next_group, 0);
            //uint64_t q = cfgMemoize (*it);
            uint64_t cost2 = cost_next_groups + cost_curr_group;
        
            //FIXME:
            /*vector <uint128_t>* _n = OptGroup::nextGroupsHashID (*it);
            if (n->size () == 0)
            {
                if (hash_id && (!(hash_id & (hash_id-1))))
                    return 0;
            }*/
            
            cost(*it, tile_sizes_next_group);
            
            if (cost1 >= 0 && cost2 >= 0)
                totalMinCost = min(totalMinCost, min (cost1, cost2));
            else if (cost1 < 0)
            {
                totalMinCost = min(totalMinCost, cost2);
            }
            else if (cost2 < 0)
            {
                totalMinCost = min(totalMinCost, cost1);
            }
            
            //if (hash_id == 1)
            if (DEBUG)
                std::cout<<"hash_id " << std::bitset<41>((uint64_t)hash_id) << " totalmin " << totalMinCost << "  " << std::bitset<41>((uint64_t)cg) << "   " << (int64_t)cost1 << "   " << cost2 << " " << std::bitset<41>((uint64_t)(*it)) << std::endl;
            if (cost1 == totalMinCost)
            {
                opt_tile_sizes = tile_sizes_cg;
                optHashID = cg;
            }
            else
            {
                opt_tile_sizes = tile_sizes_curr_group;
                opt_tile_sizes_next_group = tile_sizes_next_group;
                next_hash_id = *it;
                //optHashID = 0;
            }
        }

        if (optHashID != 0)
        {
            if (DEBUG)
                std::cout<<"hash_id " << std::bitset<41>((uint64_t)hash_id) << " opthashid " << std::bitset<41>((uint64_t)optHashID) << std::endl;
            hashIDToOptHashID[hash_id] = optHashID;
            optHashIDToTileSizes[optHashID] = opt_tile_sizes;
            
            //std::cout<< "SETTING TILE SIZES FOR " << std::bitset<41>((uint64_t)optHashID)<< std::endl;
            for (auto it = opt_tile_sizes.begin (); it != opt_tile_sizes.end (); it++)
            {
                //std::cout<< "   " << *it << std::endl;
            }
        }
        
        else{//std::cout<<"hash_id " << std::bitset<41>((uint64_t)hash_id) << " opthashid " << std::bitset<41>((uint64_t)hash_id) << std::endl;
            hashIDToOptHashID[hash_id] = hash_id;
            totalMinCost = cost_next_groups + cost_curr_group;
            optHashIDToTileSizes[hash_id] = opt_tile_sizes;
            optHashIDToTileSizes[next_hash_id] = opt_tile_sizes_next_group;
            
            /*//std::cout<< "SETTING NEXT TILE SIZES FOR " << std::bitset<41>((uint64_t)next_hash_id)<< std::endl;
            for (auto it = opt_tile_sizes_next_group.begin (); it != opt_tile_sizes_next_group.end (); it++)
            {
                //std::cout<< "   " << *it << std::endl;
            }
            
            //std::cout<< "max dimensions  " << get_max_dimensions_for_group (hash_id) << std::endl;
            uint64_t i = cost (hash_id, tile_sizes_curr_group);
            //std::cout<< "SETTING CURR TILE SIZES FOR " << std::bitset<41>((uint64_t)hash_id)<< " cost " << i << std::endl;
            for (auto it = tile_sizes_curr_group.begin (); it != tile_sizes_curr_group.end (); it++)
            {
                //std::cout<< "   " << *it << std::endl;
            }*/
        }
    }
    else
    {
        totalMinCost = 0;
        //std::cout<< hash_id << "  hash_id"<<std::endl;
        for (auto it = n->begin(); it != n->end(); ++it)
        {
            std::vector <uint64_t> tile_sizes_next_group;
            int max_dim_next_group = get_max_dimensions_for_group (*it);
            tile_sizes_next_group = std::vector <uint64_t> (max_dim_next_group, 0);
            //std::cout<<"filled it " << std::bitset<41>((uint64_t)*it)<<std::endl;
            uint64_t q = cfgMemoize (*it);
            totalMinCost += q + cost(*it, tile_sizes_next_group);
            optHashIDToTileSizes [*it] = tile_sizes_next_group;
            ////std::cout<<"filled it " << std::bitset<41>((uint64_t)*it)<<" cost "<< q + cost(*it)<<std::endl;         
        }
        
        int max_dim_curr_group = get_max_dimensions_for_group (hash_id);
        std::vector <uint64_t> tile_sizes_next_group;
        
        tile_sizes_next_group = std::vector <uint64_t> (max_dim_curr_group, 0);
        totalMinCost += cost (hash_id, tile_sizes_next_group);
        hashIDToOptHashID[hash_id] = hash_id;
        optHashIDToTileSizes[hash_id] = tile_sizes_next_group;
        //std::cout<<"hash_id " << hash_id << " cost " << totalMinCost<<std::endl;
    }
        
    T[hash_id] = totalMinCost;
    //std::cout<<"hash_id " << std::bitset <41>((uint64_t)hash_id) << " totalMin " << totalMinCost << std::endl;
    delete n;
    std::cout << "ENDING for hash_id " << std::bitset<41> ((uint64_t)hash_id) << std::endl;
    
    return totalMinCost;
}

void topological_sort (uint128_t vertex, vector<uint128_t>& output, unordered_set<uint128_t, uint128Hasher>& visited, vector<uint128_t>* next_groups)
{
    visited.insert (vertex);
    
    vector<uint128_t>* n;
    
    n = OptGroup::nextGroupsHashID (vertex);
    
    for (auto it = n->begin(); it != n->end(); it++)
    {
        if (visited.find(*it) == visited.end())
            topological_sort (*it, output, visited, next_groups);
    }
    
    if (next_groups == NULL)
        output.push_back (vertex);
        
    else if (std::find(next_groups->begin(), next_groups->end(), vertex) != next_groups->end())
    {
        output.push_back (vertex);
    }
    
    delete n;
}

bool insertingNextPrevCreatesCycle (OptGroup* g1, OptGroup* g2)
{
    bool g1_in_g2_next = g2->nextGroups().find (g1) != g2->nextGroups().end();
    bool g1_in_g2_prev = true;//g2->prevGroups().find (g1) != g2->prevGroups().end();
    bool g2_in_g1_next = g1->nextGroups().find (g2) != g1->nextGroups().end();
    bool g2_in_g1_prev = true;//g1->prevGroups().find (g2) != g1->prevGroups().end();
    
    return g1_in_g2_next && g1_in_g2_prev && g2_in_g1_next && g2_in_g1_prev;
}

OptGroup* opt_grouping_topological_sort (uint128_t _grp_hash_id)
{
    vector<uint128_t> optgroup_vector;
    unordered_set<uint128_t, uint128Hasher> visited;
    vector<uint128_t> topological_sort_output;
    //std::cout<< "SS "<<std::bitset <41>((uint64_t)_grp_hash_id)<<std::endl;
    topological_sort (_grp_hash_id, topological_sort_output, visited, NULL);
    //std::cout<<"topo done "<< std::endl;
    reverse (topological_sort_output.begin (), topological_sort_output.end ());
    
    for (auto it = topological_sort_output.begin (); it != topological_sort_output.end (); it++)
    {
        bool findInOptHashDict = false;
        //std::cout<< "going to find " << std::bitset <41>((uint64_t)*it)<<std::endl;
        for (auto it2 = optgroup_vector.begin (); it2 != optgroup_vector.end (); it2++)
        {
            if (((*it2) & (*it)) == (*it))
            {
                OptGroup* g = OptGroup::optHashIDToGroup[*it];
                //std::cout<< "exisghghts   " <<(uint64_t)g << "  " << std::bitset <41>((uint64_t)(*it))<<std::endl;
                findInOptHashDict = true;
                for (auto it3 = g->prevGroups().begin(); it3 != g->prevGroups().end(); it3++)
                {
                    for (auto it4 = (*it3)->parentGroups().begin(); it4 != (*it3)->parentGroups().end(); it4++)
                    {
                        OptGroup* optg = OptGroup::parentOptGroupsToHashID[*it2];
                        if (*it4 != optg && !insertingNextPrevCreatesCycle (*it4, optg))
                        {
                            //std::cout<<"Add " << std::bitset <41>((uint64_t)(*it4)->hashID ()) << " to " << std::bitset <41>((uint64_t)(optg)->hashID ()) << std::endl;
                            (*it4)->nextGroups().insert (optg);
                            optg->prevGroups().insert (*it4);
                        }
                    }
                }
                
                //break;
            }
        }
        
        if (findInOptHashDict == true)
            continue;
        
        uint128_t hash_id = *it;
        //std::cout<< "here    " <<std::bitset <41>((uint64_t)*it)<<std::endl;
        while (hashIDToOptHashID.find (hash_id) != hashIDToOptHashID.end())
        {
            if (hashIDToOptHashID[hash_id] == hash_id)
                break;
            
            hash_id = hashIDToOptHashID[hash_id];
        }
        //std::cout<< "SSSSSSS " <<std::bitset <41>((uint64_t)hash_id)<<std::endl;
        OptGroup* optg = OptGroup::createParentOptGroup (hash_id);
        //std::cout<< "exist in parentOptGroup " << OptGroup::parentOptGroupsToHashID [hash_id] << std::endl;
        OptGroup* g = OptGroup::optHashIDToGroup[*it];
        //std::cout<< "exists   " <<std::bitset <41>((uint64_t)(*it))<<std::endl;
        //std::cout<<"parentgroups " << (uint64_t)g << " " << g->parentGroups().size()<<std::endl;
        if (g->parentGroups().size() >= 1)
        {
            for (auto i = g->parentGroups().begin(); i != g->parentGroups().end(); i++)
            {
                //std::cout<<"parentgroup: " << std::bitset <41>((uint64_t)(*i)->hashID ()) << std::endl;
            }
        }
        
        for (auto it3 = g->prevGroups().begin(); it3 != g->prevGroups().end(); it3++)
        {
            //std::cout<<"prevGroup " << std::bitset <41>((uint64_t)(*it3)->hashID ()) << std::endl;
            for (auto it4 = (*it3)->parentGroups().begin(); it4 != (*it3)->parentGroups().end(); it4++)
            {
                if (*it4 != optg && !insertingNextPrevCreatesCycle (*it4, optg))
                {
                    //std::cout<<"Add prev" << std::bitset <41>((uint64_t)(*it4)->hashID ()) << " to " << std::bitset <41>((uint64_t)(optg)->hashID ()) << std::endl;
                    (*it4)->nextGroups().insert (optg);
                    optg->prevGroups().insert (*it4);
                }
            }
        }
        
        /*for (auto it3 = g->nextGroups().begin(); it3 != g->nextGroups().end(); it3++)
        {
            //std::cout<<"nextGroup " << std::bitset <41>((uint64_t)(*it3)->hashID ()) << std::endl;
            
            for (auto it4 = (*it3)->parentGroups().begin(); it4 != (*it3)->parentGroups().end(); it4++)
            {
                if (*it4 != optg && !insertingNextPrevCreatesCycle (*it4, optg))
                {
                    //std::cout<<"Add next " << std::bitset <41>((uint64_t)(*it4)->hashID ()) << " to " << std::bitset <41>((uint64_t)(optg)->hashID ()) << std::endl;
                    optg->nextGroups().insert (*it4);
                    (*it4)->prevGroups().insert (optg);
                }
            }
        }*/
        
        optgroup_vector.push_back (hash_id);
    }
    
    //std::cout<< "SSS DONE " << OptGroup::optHashIDToGroup[_grp_hash_id] << "        " << OptGroup::optHashIDToGroup[_grp_hash_id]->parentGroups().size() << std::endl;
    //std::cout<<"DDFDFD " << std::bitset <41>((uint64_t) (*(OptGroup::optHashIDToGroup[_grp_hash_id]->parentGroups().begin()))->hashID())<<std::endl;
    //std::cout<< "FINAL GROUP exist in parentOptGroup " << OptGroup::parentOptGroupsToHashID [(*(OptGroup::optHashIDToGroup[_grp_hash_id]->parentGroups().begin()))->hashID()] << std::endl;
        
    return *(OptGroup::optHashIDToGroup[_grp_hash_id]->parentGroups().begin());
}

OptGroup* opt_grouping (uint128_t _grp_hash_id)
{
    //std::cout<< "SSSSSSS " <<std::bitset <41>((uint64_t)_grp_hash_id)<<std::endl;
    while (hashIDToOptHashID.find (_grp_hash_id) != hashIDToOptHashID.end())
    {
        if (hashIDToOptHashID[_grp_hash_id] == _grp_hash_id)
            break;
            
        _grp_hash_id = hashIDToOptHashID[_grp_hash_id];
    }
    //std::cout<< "SS "<<std::bitset <41>((uint64_t)_grp_hash_id)<<std::endl;
    OptGroup* optGroupStart = OptGroup::createOptGroup (_grp_hash_id);
    vector<uint128_t>* next_groups = OptGroup::nextGroupsHashID (_grp_hash_id);
    vector<OptGroup*> _groups;
    unordered_set<uint128_t, uint128Hasher> visited;
    vector<uint128_t> topological_sort_output;
    //std::cout<<"topological starting "<<next_groups->size() << std::endl;
    for (auto it = next_groups->begin (); it != next_groups->end (); it++)
    {
        if (visited.find (*it) == visited.end())
            topological_sort (*it, topological_sort_output, visited, next_groups);
    }
    for (auto it = topological_sort_output.begin (); it != topological_sort_output.end (); it++)
    {
        //std::cout<< std::bitset <41>((uint64_t)_grp_hash_id) << "  " << std::bitset<41>((uint64_t)(*it)) << std::endl;
    }
    
    //std::cout<<"topological done "<<std::endl;
    reverse (topological_sort_output.begin (), topological_sort_output.end ());
    auto it = topological_sort_output.begin ();
    
    while (it != topological_sort_output.end ())
    {
        //std::cout<< "G  " <<std::bitset <41>((uint64_t)_grp_hash_id) << "  " << std::bitset<41>((uint64_t)(*it)) << std::endl;
        OptGroup* g = opt_grouping (*it);
        //std::cout<<"B   " <<std::bitset<41>((uint64_t)(g->hashID())) << "      " << std::bitset<41>((uint64_t)(*it)) << std::endl;
        bool breaked = false;
        for (; it != topological_sort_output.end (); it++)
        {
            //std::cout<<"DD      "<<std::bitset <41>((uint64_t)_grp_hash_id)<<"            "<<std::bitset<41>((uint64_t)(*it)) << std::endl;
            if ((g->hashID() & (*it)) == (*it))
            {
                //std::cout<<"CONTINUE   "<<std::bitset <41>((uint64_t)(g)->hashID())<<"        "<<std::bitset<41>((uint64_t)(*it)) << std::endl;
                continue;
            }
            else
            {
                //breaked = true;
                //std::cout<<"BREAK   "<<std::bitset <41>((uint64_t)(g)->hashID())<<"        "<<std::bitset<41>((uint64_t)(*it)) << std::endl;
                it--;
                break;
            }
        }
        
        _groups.push_back (g);
        
        if (it == topological_sort_output.end ())
            break;
        //if (breaked == false)
            it++;
    }
    
    /*size_t prev_size = 0;
    
    while (prev_size != _groups.size ())
    {
        prev_size = _groups.size ();
        
        vector<OptGroup*>::iterator to_remove = _groups.end();
        
        for (auto it = _groups.begin(); it != _groups.end (); it++)
        {
            for (auto it2 = _groups.begin (); it2 != _groups.end (); it2++)
            {
                if ((*it)->hashID() != (*it2)->hashID() && (*it)->hashID() & (*it2)->hashID() == (*it)->hashID())
                {
                    to_remove = it;
                    break;
                }
            }
            
            if (to_remove != _groups.end())
                break;
        }
        
        if (to_remove != _groups.end())
        {
            //std::cout<<"prev_size "<<prev_size << " to remove "<< std::bitset<41>((uint64_t)(*to_remove)->hashID())<<std::endl;
            _groups.erase (to_remove);
            OptGroup::Q.erase ((*to_remove)->hashID());
            (*to_remove)->remove_from_parent ();
        }
    }*/
    
    for (auto it = _groups.begin(); it != _groups.end(); it++)
    {
        optGroupStart->nextGroups().insert(*it);
    }
    
    vector<OptGroup*> to_remove;
    //std::cout<<"determine to_remove"<<std::endl;
    for (auto it2 = optGroupStart->nextGroups().begin (); it2 != optGroupStart->nextGroups().end (); it2++)
    {
        uint128_t hash_id = (*it2)->hashID();
        
        if (hash_id && (!(hash_id & (hash_id-1))))
        {
            for (auto it = optGroupStart->nextGroups().begin (); it != optGroupStart->nextGroups().end (); it++)
            {
                if (hash_id && (*it)->hashID() == hash_id)
                {
                    to_remove.push_back ((OptGroup*)*it2);
                }
            }
        }
    }
    
    //std::cout<<"to_remove vertices"<<std::endl;
    
    for (auto it = to_remove.begin(); it != to_remove.end(); it++)
    {
        optGroupStart->nextGroups().erase (*it);
    }
    
    //std::cout<<"returning optGroupStart"<<std::endl;
    
    return optGroupStart;
}

bool isCyclicUtil (OptGroup* v, unordered_map<OptGroup*, bool>& visited, 
                   unordered_map<OptGroup*, bool>& recStack, 
                   OptGroup** cycle, int cycleIndex, int& cycle_length)
{
    //std::cout<< "isCyclicUtil called " << std::endl;
    if(visited[v] == false)
    {
        //std::cout<< "isCyclicUtil " << cycle_length << std::endl;
        // Mark the current node as visited and part of recursion stack
        visited[v] = true;
        recStack[v] = true;
        cycle [cycleIndex] = v;
        // Recur for all the vertices adjacent to this vertex
        for(auto &iter : v->nextGroups ())
        {
            cycle_length ++;
            if (!visited[(OptGroup*)iter] && isCyclicUtil((OptGroup*)iter, visited, recStack, 
                                               cycle, cycleIndex + 1, 
                                               cycle_length))
            {
                return true;
            }
            else if (recStack[(OptGroup*)iter])
            {
                return true;
            }
            cycle_length--;
        }
    }
    
    recStack[v] = false;  // remove the vertex from recursion stack
    return false;
}

void detect_and_destroy_cycles ()
{
    bool change;
    
    change = true;
    
    while (change)
    {
        // Mark all the vertices as not visited and not part of recursion
        // stack
        int n_vertices;
        unordered_map <OptGroup*, bool> visited;
        unordered_map <OptGroup*, bool> recStack;
        n_vertices = OptGroup::parentOptGroupsToHashID.size ();
        OptGroup* cycle [n_vertices];
        int cycle_length = 0;
        
        change = false;
        memset (cycle, sizeof (OptGroup*)*n_vertices, 0);
        
        for(auto const &iter: OptGroup::parentOptGroupsToHashID)
        {
            visited[iter.second] = false;
            recStack[iter.second] = false;
        }
        static int n_calls = 0;
        n_calls++;
        //std::cout<< "change is " << change << std::endl;
        // Call the recursive helper function to detect cycle in different
        // DFS trees
        int q = 0;
        for(auto const &iter: OptGroup::parentOptGroupsToHashID)
        {
            q++;
            cycle_length = 0;
            
            if (isCyclicUtil(iter.second, visited, recStack, cycle, 0, cycle_length))
            {
                //std::cout<<"Found Cycle #" <<q << " in call #"<< n_calls << std::endl;
                for (int k = 0; k < cycle_length; k++)
                    //std::cout<< std::bitset<41>((uint64_t)cycle [k]->hashID ()) << std::endl;
                
                change = true;
                break;
            }
        }
        
        if (!change)
            continue;
            
        //Break cycle by creating separate group of child which is not the input (of group)
        //and has an edge from another group
        OptGroup* child_to_disconnect = NULL;
        OptGroup* curr_grp;
        OptGroup* next_grp;
        
        for (int i = 0; i < cycle_length; i++)
        {
            int curr_grp_index = i;
            int next_grp_index = (i+1)%cycle_length;
            curr_grp = cycle [curr_grp_index];
            next_grp = cycle [next_grp_index];
            
            //Find a child in next_grp having an edge from a child of curr_grp
            //and is not the input of the next_grp
            for (auto const &curr_child: curr_grp->childGroups ())
            {
                for (auto const &next_curr_child: curr_child->nextGroups ())
                {
                    if ((next_curr_child->hashID () & next_grp->hashID ()) ==  
                        (next_curr_child->hashID ()))
                    {
                       bool candidate_is_input = false;
                       
                       for (auto const &prev : next_curr_child->prevGroups ())
                       {
                           if ((prev->hashID () & next_grp->hashID ()) ==
                                prev->hashID ())
                           {
                               candidate_is_input = true;
                               break;
                           }
                       }
                       
                       //std::cout<< "candidate_is_input " << candidate_is_input << std::endl;
                       if (candidate_is_input)
                            child_to_disconnect = (OptGroup*)next_curr_child;
                    }
                
                    if (child_to_disconnect)
                        break;
                }
                
                if (child_to_disconnect)
                    break;
            }
            
            if (child_to_disconnect)
                break;
        }
        
        //std::cout<< "child_to_disconnect = " <<child_to_disconnect << std::endl;
        //std::cout<< "child_to_disconnect is " << std::bitset<41>((uint64_t)child_to_disconnect->hashID()) << std::endl;
        
        OptGroup* to_remove = OptGroup::createParentOptGroup (child_to_disconnect->hashID ());
        //std::cout<< "to_remove " << to_remove << std::endl;
        
        assert (to_remove != child_to_disconnect);
        
        for (auto const &_next_grp: child_to_disconnect->nextGroups ())
        {
            for (auto const &grp : ((OptGroup*)_next_grp)->parentGroups ())
            {
                if (grp->hashID () != next_grp->hashID ())
                    to_remove->nextGroups ().insert (grp);
            }
        }
        
        for (auto const &_next_grp : to_remove->nextGroups ())
        {
            if (_next_grp->hashID () != next_grp->hashID ())
                _next_grp->prevGroups ().insert ((OptGroup*)to_remove);
        }
        
        //std::cout<< "next_grp_done" << std::endl;
        for (auto const &prev_grp: child_to_disconnect->prevGroups ())
        {
            for (auto const &grp : ((OptGroup*)prev_grp)->parentGroups ())
            {
                if (grp->hashID () != next_grp->hashID ())
                    to_remove->prevGroups ().insert (grp);
            }
        }
        
        for (auto const &prev_grp : to_remove->prevGroups ())
        {
            if (prev_grp->hashID () != next_grp->hashID ())
                prev_grp->nextGroups ().insert ((OptGroup*)to_remove);
        }
        
        //std::cout<< "prev_grp_done " << std::endl;
        
        //Insert to_remove in prevGroups and nextGroups of 
        //to_remove->nextGroups and to_remove->prevGroups
        next_grp->removeChild (to_remove->hashID (), to_remove);
        //std::cout<< "232323 " << std::endl;
    }
}

void correct_grouping ()
{
    bool change = true;
    
    while (change)
    {
        change = false;
        for (auto it = OptGroup::optHashIDToGroup.begin (); 
             it != OptGroup::optHashIDToGroup.end (); it++)
        {
            OptGroup* group = it->second;
            
            if (group == NULL)
                continue;
            
            if (DEBUG){    
            std::cout<<"LLLLL "<< std::bitset<41>((uint64_t)group->hashID ())<< std::endl;
            std::cout<<"N_PARENTS " << group->parentGroups().size() << std::endl;
            }
        
            for (auto it2 = group->parentGroups().begin(); 
                 it2 != group->parentGroups().end(); it2++)
            {
                if (DEBUG){
                    std::cout<<"PARENT " << std::bitset<41>((uint64_t)(*it2)->hashID())<<std::endl;
                }
            }
            
            if (group != NULL && group->parentGroups().size () <= 1)
                continue;
    
            change = true;
            bool in_middle = false;
            OptGroup* in_middle_group = 0;
            
            //sort (group->parentGroups().begin (), group->parentGroups().end(), optGroupVecFunc);
            
            /*for (auto it = group->parentGroups().begin (); it != group->parentGroups().end(); it++)
            {
                uint64_t first_1_pos = first_one ((*it)->hashID());
                uint64_t n_ones = numberOfOnes ((*it)->hashID()) - 1;
                if (((uint64_t)id) != n_ones + first_1_pos && ((uint64_t)id) != first_1_pos)
                {
                    in_middle = true;
                    in_middle_group = *it;
                    break;
                }
            }
            
            for (auto it = group->parentGroups().begin (); it != group->parentGroups().end(); it++)
            {
                if (*it != in_middle_group)
                {
                    (*it)->removeChild (group->hashID());
                }
            }
            
            if (in_middle)
                continue;
            */
            OptGroup* minParentGroup = 0;
            uint64_t minCost = 1L << 60;
            
            for (auto it = group->parentGroups().begin (); 
                 it != group->parentGroups().end(); it++)
            {
                //std::cout<<"parent groups: "<< std::bitset<41>((uint64_t)(*it)->hashID ())<< std::endl;
                int max_dim = get_max_dimensions_for_group ((*it)->hashID());
                uint64_t _cost;
                std::vector <uint64_t> tile_sizes;
                
                tile_sizes = std::vector <uint64_t> (max_dim, 0); 
                _cost = cost ((*it)->hashID(), tile_sizes);
                
                if (_cost < minCost && !((*it)->hashID() && (!((*it)->hashID() & ((*it)->hashID()-1)))))
                {
                    minCost = _cost;
                    minParentGroup = (*it);
                }
            }
            
            if (DEBUG){
            std::cout<<"wewwewe "<< std::bitset<41>((uint64_t)group->hashID ())<< std::endl;
            std::cout<<"minParentGroup " << std::bitset<41>((uint64_t)minParentGroup->hashID())<< std::endl;
            }
        
            /* removeChild () modifies the child's parentGroups() set.
             * Hence, we cannot call removeChild () from the loop iterator 
             * of child's parentGroups ().
             * Create a copy of parents and call removeChild on them.
             * */
            vector<OptGroup*> toRemoveParents;
            
            for (auto it = group->parentGroups().begin (); 
                 it != group->parentGroups().end(); it++)
            {
                if (*it != minParentGroup)
                {
                    toRemoveParents.push_back (*it);
                }
            }
            
            for (auto it = toRemoveParents.begin(); it != toRemoveParents.end();
                 it++)
            {   
                if (DEBUG){
                    std::cout<<"removechild "<< std::endl;
                }
                (*it)->removeChild (group->hashID(), minParentGroup);
            }
            if (DEBUG){
                std::cout<<"wew2222 "<< std::bitset<41>((uint64_t)group->hashID ())<< std::endl;
            }
            
            /*uint128_t a = 0b00001000000000000000000000000000000000000;
            OptGroup* gggg = OptGroup::parentOptGroupsToHashID[a];
            if (gggg)
            {
                //std::cout<<"prevGroups for buggy group " << std::endl;
                for (auto const &g : gggg->prevGroups ())
                {
                    //std::cout<< "   " << std::bitset <41>((uint64_t)(g->hashID ())) << std::endl;
                }
                int ___i;
                std::cin >> ___i;
                std::cin >> ___i;
            }*/
        }
    }
    
    //std::cout<< "corrected grouping before detecting cycles"<<std::endl;
        for (auto it = OptGroup::parentOptGroupsToHashID.begin(); it != OptGroup::parentOptGroupsToHashID.end(); it++)
        {   
            if (it->second == NULL)
                continue;
                
            //std::cout<< std::bitset <41>((uint64_t)it->second->hashID ()) << std::endl;
            ////std::cout<< it->second->hashID() << std::endl;
            for (auto it2 = it->second->nextGroups().begin(); it2 != it->second->nextGroups().end(); it2++)
            {
                //std::cout<< "        "<< std::bitset <41>((uint64_t)(*it2)->hashID ()) << std::endl;
                ////std::cout<< "        "<< (*it2)->hashID () << std::endl;
            }
        }
    
    //std::cout<< "corrected grouping prevGroups before detecting cycles"<<std::endl;
        for (auto it = OptGroup::parentOptGroupsToHashID.begin(); it != OptGroup::parentOptGroupsToHashID.end(); it++)
        {   
            if (it->second == NULL)
                continue;
                
            //std::cout<< std::bitset <41>((uint64_t)it->second->hashID ()) << std::endl;
            ////std::cout<< it->second->hashID() << std::endl;
            for (auto it2 = it->second->prevGroups().begin(); it2 != it->second->prevGroups().end(); it2++)
            {
                //std::cout<< "        "<< std::bitset <41>((uint64_t)(*it2)->hashID ()) << std::endl;
                ////std::cout<< "        "<< (*it2)->hashID () << std::endl;
            }
        }
    //FIXME: Cycles are created with __SERVER__ configuration with cores = 16
    //Below is just a temporary fix, to keep the project progressing 
    detect_and_destroy_cycles ();
}

PyObject* firstPyGroup (uint128_t hash_id, bool dummySource, uint64_t maxID)
{
    uint64_t bit = 0;
        
    while (hash_id != 0)
    {
        if (hash_id & 1L == 1)
        {   
            uint128_t l = 1;
            l = l<<bit;
            //std::cout<<" firstpygroup hash id "<< std::bitset<41>((uint64_t)l) << std::endl;
            Group* g = Group::hashIDToGroup[l];
            //std::cout<<"firstpygroup "<< g << std::endl;
            PyObject* o = g->getPyGroup ();
            if (o != NULL)
                return o;
        }
        
        hash_id = hash_id >> 1;
        bit++;
    }
    //assert (false);
    return NULL;

}

PyObject* getPyGroupTopologicalOrder (PyObject* pyg)
{
    PyObject* order = PyDict_GetItem (pygroup_topological_order, 
                                      pyg);
    return order;
}

bool pyGroupSortFunc (PyObject* pyg1, PyObject* pyg2)
{    
    return PyLong_AsLong (getPyGroupTopologicalOrder (pyg1)) < PyLong_AsLong (getPyGroupTopologicalOrder (pyg2));
}

void clear_everything ()
{
    hashIDToOptHashID.clear();
    T.clear();
    vector<std::unordered_map <uint128_t, OptGroup*, uint128Hasher>::iterator> toremove;
    
    OptGroup::optHashIDToGroup.clear ();
    OptGroup::vectorOptGroups.clear();
    
    for (auto it = OptGroup::parentOptGroupsToHashID.begin(); it!= OptGroup::parentOptGroupsToHashID.end(); it++)
    {
        OptGroup::optHashIDToGroup[it->first] = it->second;
        OptGroup::vectorOptGroups.push_back (it->second);
    }
    
    OptGroup::parentOptGroupsToHashID.clear ();
    
    optHashIDToTileSizes.clear ();
    
}

/*
 * Main function of the library. Arguments are 
 * in_group: A list of input groups
 * out_group: A list of output groups
 * groups: A list of all groups
 * pipeline: The Pipeline Object
 * reduction_cls: Class object of Reduction class
 * small_comps: A list of small computations
 * comp_size_map: A map of each computation to size
 * tstencil_cls: Class object of TStencil Class
 * pygroup_topological_order: Topological Order of Groups
 * pygroup_dim_reuse: Dimension reuse across each dimension for each group
 * pylive_size: 
 * pydim_size: Dimension sizes of each dimension of each group
 * storage_mapping_get_dim_size: Object for function storage_mapping.get_dim_size
 * cls_Storage: Class Object of Storage Class
 */
PyObject* py_group(PyObject* self, PyObject* args)
{
    bool dummySource;
    PyObject *in_group;
    PyObject *out_group;
    PyObject *groups;    
    
    std::vector <PyObject*> pygroups_vector;
    
    std::cout<<"PolyMage+ grouping " << std::endl;
    int ee;
    std::cin >> ee;
    std::cin >> ee;
    PyArg_ParseTuple (args, "OOOOOOOOOOOOOO", &in_group, &out_group, &groups, &pipeline, 
                      &reduction_cls, &small_comps, &comp_size_map, &tstencil_cls,
                      &pygroup_topological_order, &pygroup_dim_reuse, &pylive_size,
                      &pydim_size, &storage_mapping_get_dim_size, &cls_Storage);
    get_overlapping_size_func = PyObject_GetAttr (pipeline, 
                                                  Py_BuildValue ("s", 
                                                  "get_overlapping_size_for_groups"));
    
    if (PyList_Size (in_group) == 1)
    {
        maxID = PyList_Size (groups);
        dummySource = false;
    }
    else
    {
        maxID = PyList_Size (groups) + 1;
        dummySource = true;
    }

    if (DEBUG)
    {
        std::cout<< "Total Number of Nodes "<< maxID << std::endl;
    }
    
    for (int i = 0; i < PyList_Size (small_comps); i++)
    {
        PyObject* comp = PyList_GetItem (small_comps, i);
        
        small_comps_set.insert (comp);
    }
    
    for (int i = 0; i < PyList_Size (groups); i++)
    {
        pygroups_vector.push_back (PyList_GetItem (groups, i));
    }
    
    sort (pygroups_vector.begin (), pygroups_vector.end(), pyGroupSortFunc);
    
    for (int i = 0; i < pygroups_vector.size(); i++)
    {
        PyObject* group = pygroups_vector [i];
        PyObject* comps = PyObject_GetAttr (group,
                                            Py_BuildValue ("s", 
                                                           "comps"));
        for (int j = 0; j < PyList_Size (comps); j++)
        {
            static PyObject* str_lookup_key = Py_BuildValue ("s", "lookup_key");
            static PyObject* str_has_key = Py_BuildValue ("s", "has_key");
            static PyObject* str_orig_stg_class = Py_BuildValue ("s", "orig_storage_class");
            static PyObject* str_func = Py_BuildValue ("s", "func");
            static PyObject* str_typ = Py_BuildValue ("s", "typ");
            static PyObject* str_ndims = Py_BuildValue ("s", "ndims");
            static PyObject* str_offsets = Py_BuildValue ("s", "offsets");
            static PyObject* str_generate_id = Py_BuildValue ("s", "generate_id");
            PyObject* storage;
            PyObject* key;
            PyObject* comp;
            PyObject* func;
            PyObject* offsets;
            PyObject* typ;
            PyObject* dims;
            
            comp = PyList_GetItem (comps, j);
            func = PyObject_GetAttr (comp, str_func);
            typ = PyObject_GetAttr (func, str_typ);
            dims = PyObject_GetAttr (func, str_ndims);
            storage = PyObject_GetAttr (comp, str_orig_stg_class);
            key = PyObject_GetAttr (storage, str_lookup_key);
            string _strkey = pystring_to_string (key);
        
            pyGroupForComp [comp] = group;
            comp_to_orig_storage [comp] = storage;
            stg_class_to_key [storage] = _strkey;
            comp_to_func [comp] = func;
            func_to_typ [func] = typ;
            func_to_ndims [func] = PyLong_AsLong (dims);
            
            offsets = PyObject_GetAttr (storage, str_offsets);
            //check_and_print_exception ();
            vector <long> max_offset;
            
            for (long i = 0; i < PyLong_AsLong (dims); i++)
            {
                PyObject* _max_offset = PyTuple_GetItem (PyList_GetItem (offsets, i), 1);
                //check_and_print_exception ();
                max_offset.push_back (PyLong_AsLong (_max_offset));
                //check_and_print_exception ();
            }
            
            comp_to_max_offsets [comp] = max_offset;
            helper_storage_to_dim_storage [comp] = std::vector<PyObject*> ();
            
            for (long dim = 0; dim < PyLong_AsLong (dims); dim++)
            {
                static PyObject* str_get_dim = Py_BuildValue ("s", "get_dim");
                static PyObject* str_coeff = Py_BuildValue ("s", "coeff");
                static PyObject* str_orig_param = Py_BuildValue ("s", "orig_param");
                PyObject* dim_storage;
                PyObject* get_dim;
                PyObject* new_size;
                
                get_dim = PyObject_GetAttr (storage, 
                                            str_get_dim);
                //check_and_print_exception ();
                PyObject* args = PyTuple_Pack (1, PyLong_FromLong (dim));
                //check_and_print_exception ();
                dim_storage = PyObject_CallObject (get_dim, args);
                //check_and_print_exception ();
                Py_DECREF (args);
                
                helper_storage_to_dim_storage [storage].push_back (dim_storage);
            }
        }
        
        Group* cppGroup = new Group (group, comps);
        pyToGroup[group] = cppGroup;
        groupToPyGroup[cppGroup] = group;
        OptGroup* _opt_group = OptGroup::createOptGroup (cppGroup->hashID ());
        PyObject* name = PyObject_GetAttr (group,
                                                   Py_BuildValue ("s", 
                                                                  "name"));
        PyObject* pyStr = PyUnicode_AsEncodedString(name, "utf-8", "Error ~");
        const char* ss = PyBytes_AS_STRING (pyStr);
        
        if (DEBUG)
        {
            std::cout<< "hash_id " << std::bitset<41>((uint64_t)cppGroup->hashID()) << "  " << ss << std::endl;
        }
    }
    
    if (DEBUG)
    {
        std::cout << "comp_to_stg_class size " << comp_to_orig_storage.size() << std::endl;
        for (auto it = comp_to_orig_storage.begin(); it != comp_to_orig_storage.end(); it++)
        {
            std::cout << it->first << "  " << it->second <<std::endl;
        }
    }
    
    for (int i = 0; i < PyList_Size (groups); i++)
    {
        PyObject* group = PyList_GetItem (groups, i);
        Group* cppGroup = pyToGroup[group];
        OptGroup* optGroup = OptGroup::optHashIDToGroup [cppGroup->hashID ()];
        PyObject* children = PyObject_GetAttr (group,
                                               Py_BuildValue ("s", 
                                                              "children"));        
        for (int j = 0; j < PyList_Size (children); j++)
        {
            PyObject* child = PyList_GetItem (children, j);
            Group* cppChild = pyToGroup [child];
            cppGroup->nextGroups().insert (cppChild);
            OptGroup* optChild = OptGroup::optHashIDToGroup [cppChild->hashID ()];
            optGroup->nextGroups().insert (optChild);
        }
        
        PyObject* parents = PyObject_GetAttr (group,
                                               Py_BuildValue ("s", 
                                                              "parents"));
        
        for (int j = 0; j < PyList_Size (parents); j++)
        {
            PyObject* parent = PyList_GetItem (parents, j);
            Group* cppParent = pyToGroup [parent];
            cppGroup->prevGroups().insert (cppParent);
            OptGroup* optParent = OptGroup::optHashIDToGroup [cppParent->hashID ()];
            optGroup->prevGroups().insert (optParent);
        }
    }
    
    Group* start;
    OptGroup* opt_start;
    
    if (PyList_Size (in_group) == 1)
    {
        start = pyToGroup[PyList_GetItem (in_group, 0)];
        opt_start = OptGroup::optHashIDToGroup[start->hashID()];
    }
    else
    {
        start = new Group ();
        opt_start = OptGroup::createOptGroup (start->hashID ());
        for (int i = 0; i < PyList_Size (in_group); i++)
        {
            Group* next = pyToGroup[PyList_GetItem (in_group, i)];
            start->nextGroups().insert (next);
            next->prevGroups().insert (start);
            OptGroup* optnext = OptGroup::optHashIDToGroup [next->hashID ()];
            opt_start->nextGroups().insert (optnext);
            optnext->prevGroups().insert (opt_start);
        }
    }
    
    PyObject* grp_size_o = PyObject_GetAttr (pipeline,
                                             Py_BuildValue ("s", 
                                                            "_group_size"));
    grp_size = PyLong_AsLong (grp_size_o);
    
    if (DEBUG)
    {
        std::cout << "Group Size " << grp_size << std::endl;
    }
    
    int max_group_size = 1 << logMaxChildren;
    uint128_t start_hash_id = opt_start->hashID();
       
    int iteration = 0;
    uint64_t _maxID = maxID;
    if (_maxID && (!(_maxID & (_maxID-1))))
    {}
    else
    {
        //Next power of 2 of maxID
        _maxID--;
        _maxID |= _maxID >> 1;
        _maxID |= _maxID >> 2;
        _maxID |= _maxID >> 4;
        _maxID |= _maxID >> 8;
        _maxID |= _maxID >> 16;
        _maxID++;
        if (DEBUG)
        {
            std::cout<<"New Number of Nodes" << _maxID<<std::endl; 
        }
    }
    
    while (iteration == 0 || max_group_size <= _maxID)
    {
        uint128_t k = 1 << 16;
        ////std::cout<< "has parents group " << (uint64_t) OptGroup::parentOptGroupsToHashID [k] << std::endl;
        //std::cout<<"clearing everythinhg"<<std::endl;
        if (iteration != 0)
            clear_everything ();
        k = 1 << 16;
        //std::cout<< "has missing group " << (uint64_t) OptGroup::optHashIDToGroup [k] << std::endl;
    
        //std::cout<< "starting cfg memoize"<<std::bitset <41>((uint64_t)start_hash_id) <<std::endl;
        uint64_t minCost = cfgMemoize (start_hash_id);
        //std::cout<< "cfg memoize ended"<<std::endl;
        
        //std::cout<< "hashIDToOptHashID"<<std::endl;
        for (auto it = hashIDToOptHashID.begin(); it != hashIDToOptHashID.end(); it++)
        {
            //std::cout<< std::bitset <41>((uint64_t)it->first) << "    " << std::bitset <41>((uint64_t)it->second)<<std::endl;
        }
        //std::cout<< "hashIDToOptHashID22222"<<std::endl<<std::endl;
        //std::cin >> iteration;
        //opt_start = opt_grouping (start->hashID());
        opt_start = opt_grouping_topological_sort (start_hash_id);
        if (DEBUG)
            std::cout<< "opt_grouping ended iteration "<<iteration<<std::endl;
        ////std::cout<<"Group " << std::endl;
        for (auto it = OptGroup::parentOptGroupsToHashID.begin(); it != OptGroup::parentOptGroupsToHashID.end(); it++)
        {
            if (it->second == NULL)
                continue;
                
            if (DEBUG)
            std::cout<< std::bitset <41>((uint64_t)it->second->hashID()) << std::endl;
            for (auto it2 = it->second->nextGroups().begin(); it2 != it->second->nextGroups().end(); it2++)
            {
                if (DEBUG)
                    std::cout<< "        "<< std::bitset <41>((uint64_t)(*it2)->hashID()) << std::endl;
            }
        }
    
        if (DEBUG)
            std::cout<< "opt_grouping prev groups"<<std::endl;
        ////std::cout<<"Group " << std::endl;
        for (auto it = OptGroup::parentOptGroupsToHashID.begin(); it != OptGroup::parentOptGroupsToHashID.end(); it++)
        {
            if (it->second == NULL)
                continue;
                
            if (DEBUG)
                std::cout<< std::bitset <41>((uint64_t)it->second->hashID()) << std::endl;
            for (auto it2 = it->second->prevGroups().begin(); it2 != it->second->prevGroups().end(); it2++)
            {
                if (DEBUG)
                    std::cout<< "        "<< std::bitset <41>((uint64_t)(*it2)->hashID()) << std::endl;
            }
        }
        
        correct_grouping();
    
        if (DEBUG)
            std::cout<< "corrected grouping"<<std::endl;
        for (auto it = OptGroup::parentOptGroupsToHashID.begin(); it != OptGroup::parentOptGroupsToHashID.end(); it++)
        {   
            if (it->second == NULL)
                continue;
                
            if (DEBUG)
                std::cout<< std::bitset <41>((uint64_t)it->second->hashID ()) << std::endl;
            ////std::cout<< it->second->hashID() << std::endl;
            for (auto it2 = it->second->nextGroups().begin(); it2 != it->second->nextGroups().end(); it2++)
            {
                if (DEBUG)
                    std::cout<< "        "<< std::bitset <41>((uint64_t)(*it2)->hashID ()) << std::endl;
                ////std::cout<< "        "<< (*it2)->hashID () << std::endl;
            }
        }
        
        if (DEBUG)
            std::cout<< "corrected grouping prevGrous"<<std::endl;
        for (auto it = OptGroup::parentOptGroupsToHashID.begin(); it != OptGroup::parentOptGroupsToHashID.end(); it++)
        {   
            if (it->second == NULL)
                continue;
                
            if (DEBUG)
                std::cout<< std::bitset <41>((uint64_t)it->second->hashID()) << std::endl;
            for (auto it2 = it->second->prevGroups().begin(); it2 != it->second->prevGroups().end(); it2++)
            {
                if (DEBUG)
                    std::cout<< "        "<< std::bitset <41>((uint64_t)(*it2)->hashID()) << std::endl;
            }
        }
        
        start_hash_id = opt_start->hashID();
        iteration++;
        logMaxChildren += log_increment_factor;
        max_group_size = 1 << logMaxChildren;
    }
    
    //std::cout<< "TOTAL SIZE USED: " << std::endl;
    for (auto it = OptGroup::parentOptGroupsToHashID.begin (); 
         it != OptGroup::parentOptGroupsToHashID.end (); it++)
    {
        //std::cout<< std::bitset <41>((uint64_t)it->second->hashID ()) << "  " << getTotalSizeUsed (it->second->hashID ()) << "Bytes" << std::endl;
    }
    
    /*//std::cout<< "TILE SIZES ARE: " << std::endl;
    for (auto it = optHashIDToTileSizes.begin (); it != optHashIDToTileSizes.end (); it++)
    {
        //std::cout<< std::bitset <41>((uint64_t)it->first) << std::endl;
        for (auto it2 = it->second.begin (); it2 != it->second.end (); it2++)
        {
            //std::cout<< "  " << *it2 << std::endl;
        }
    }*/

    
    std::cin>>iteration;
    //std::cout<<"starting bfs"<<std::endl;
    queue<OptGroup*> opt_queue;
    OptGroup* opt_temp = opt_start;
    opt_queue.push (opt_temp);
    PyObject* merge_groups_func = PyObject_GetAttr (pipeline,
                                              Py_BuildValue ("s", 
                                                             "merge_groups"));
    
    PyObject* set_tile_size_func;
                                                             
    //std::cout<<"bfs running"<<std::endl;
    std::unordered_map<uint128_t, bool, uint128Hasher> visited;
    visited[opt_temp->hashID()] = true;
    std::unordered_set<OptGroup*> done;
    //std::cout<<"opt_temp " << std::bitset <41>((uint64_t)opt_temp->hashID()) << std::endl;
    while (opt_queue.size () != 0)
    {
        opt_temp = opt_queue.front ();
        opt_queue.pop();
        /*//std::cout<< "TILE SIZEs are " << std::endl;
        for (auto it = optHashIDToTileSizes [opt_temp->hashID ()].begin (); 
             it != optHashIDToTileSizes [opt_temp->hashID ()].end (); it++)
        {
            //std::cout<< *it << std::endl;
        }*/
        
        PyObject* new_grp = firstPyGroup (opt_temp->hashID(), dummySource, maxID);
        //std::cout<<"firstPyGroup " << new_grp << "   " << std::bitset <41>((uint64_t)opt_temp->hashID()) << std::endl;
        
        if (done.find(opt_temp) == done.end() && opt_temp->hashID() != 0)
        {
            done.insert(opt_temp);
            uint64_t first = first_one (opt_temp->hashID());
            if (numberOfOnes (opt_temp->hashID ()) > 1)
            {
                //std::cout<< "sfsdf3333" << std::endl;
                uint64_t bit = 0;
                uint128_t hash_id = opt_temp->hashID();
                
                while (bit <= first)
                {
                    hash_id = hash_id >> 1;
                    bit++;
                }
                
                while (hash_id != 0)
                {                    
                    if ((hash_id & 1L) == 1 && !(dummySource && bit == maxID-1))
                    {   
                        uint128_t l = 1;
                        l = l<<bit;
                        
                        //std::cout<< "bit " << bit <<" child "<<std::bitset <41>((uint64_t)l)<<std::endl;

                        Group* g = Group::hashIDToGroup[l];
                        //std::cout<<"sfsdfsdf "<<std::bitset <41>((uint64_t)hash_id)<<std::endl;
                        PyObject* pp = g->getPyGroup ();
                        if (pp == NULL)
                            continue;
                        //std::cout<<"2323223 "<< pp <<std::endl;
                        PyObject* args = PyTuple_Pack (2, new_grp, pp);
                        ////std::cout<<"123123123123123"<<std::endl;
                        new_grp = PyObject_CallObject (merge_groups_func, args);
                        if (PyErr_Occurred () != NULL)
                        {
                            ////std::cout<<"Exception occurred "<<std::endl;
    
                            PyErr_Print();
                        }
                        //else
                            ////std::cout<<"NO Exception"<<std::endl;
                        ////std::cout<<"wwwwwwwwwwwwwwwwww "<<new_grp << std::endl;
                    }
            
                    hash_id = hash_id >> 1;
                    bit++;
                }
            }
        }
        
        ////std::cout<< "new_grp is " << new_grp <<std::endl;
        if (new_grp)
        {
            uint64_t total_size_used;
            static PyObject* set_total_used_size = Py_BuildValue ("s", 
                                                                  "set_total_used_size");
            static PyObject* set_n_buffers = Py_BuildValue ("s", 
                                                            "set_n_buffers");

            PyObject* set_total_used_size_func, *set_n_buffers_func;
            int n_buffers=0, liveout_buffers=0;
            uint64_t liveouts_size, livein_size;
            total_size_used = getTotalSizeUsed (opt_temp->hashID(), n_buffers, liveouts_size, livein_size);
            set_total_used_size_func = PyObject_GetAttr (new_grp,
                                                         set_total_used_size);
            PyObject* args = PyTuple_Pack (1,
                                           PyLong_FromLong (total_size_used));
            PyObject_CallObject (set_total_used_size_func, args);
            Py_DECREF (args);
            
            set_n_buffers_func = PyObject_GetAttr (new_grp,
                                                   set_n_buffers);
            args = PyTuple_Pack (1,
                                 PyLong_FromLong (n_buffers));
            PyObject_CallObject (set_n_buffers_func, args);
            Py_DECREF (args);
        }
        
        for (auto it = opt_temp->nextGroups().begin(); 
             it != opt_temp->nextGroups().end(); it++)
        {
            if (visited.find(((*it)->hashID())) == visited.end())
                visited[(*it)->hashID()] = false;
                
            if (visited[(*it)->hashID()] == false)
            {
                visited[opt_temp->hashID()] = true;
                opt_queue.push ((OptGroup*)*it);
            }
        }
    }
    
    pyGroupForComp.clear ();
    pyToGroup.clear ();
    groupToPyGroup.clear ();
    OptGroup::parentOptGroupsToHashID.clear ();
    hashIDToOptHashID.clear ();
    Group::hashIDToGroup.clear ();
    T.clear ();
    OptGroup::optHashIDToGroup.clear ();
    small_comps_set.clear ();
    optHashIDToTileSizes.clear ();
    pyGroupForComp.clear ();
    globalID = -1;
    //std::cout<<"bfs finished"<<std::endl;
    return Py_BuildValue ("s", "_group_size");
}

PyMethodDef myModule_methods[] = {
  {"group", py_group, METH_VARARGS},
  {NULL, NULL}
};

struct PyModuleDef myModule =
{
    PyModuleDef_HEAD_INIT,
    "optgrouping", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    myModule_methods
};

/*
 * Python calls this to let us initialize our module
 */
PyMODINIT_FUNC PyInit_optgrouping()
{
  return PyModule_Create(&myModule);
}
