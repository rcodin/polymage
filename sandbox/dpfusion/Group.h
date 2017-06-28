#include <iostream>
#include <Python.h>
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

using boost::multiprecision::uint128_t;
using namespace boost;
using namespace std;

#ifndef __GROUP_H__
#define __GROUP_H__

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
    
    static void liveoutsForGroup (const uint128_t hash_id, 
                                  const vector<Group*>& children, 
                                  unordered_set<Group*, Group::GroupHasher>& liveouts)
    {
        for (auto const& child: children)
        {
            bool is_liveout = false;
            
            for (auto const& next: child->nextGroups ())
            {
                if ((next->hashID() & hash_id) != next->hashID())
                {
                    is_liveout = true;
                    break;
                }                
            }
            
            if (child->nextGroups ().size () == 0)
                is_liveout = true;
                
            if (is_liveout)
            {
                liveouts.insert (child);
            }
        }
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

#endif
