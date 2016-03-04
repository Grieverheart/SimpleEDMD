/*
* Copyright (c) 2009 Erin Catto http://www.box2d.org
*
* This software is provided 'as-is', without any express or implied
* warranty.  In no event will the authors be held liable for any damages
* arising from the use of this software.
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely, subject to the following restrictions:
* 1. The origin of this software must not be misrepresented; you must not
* claim that you wrote the original software. If you use this software
* in a product, an acknowledgment in the product documentation would be
* appreciated but is not required.
* 2. Altered source versions must be plainly marked as such, and must not be
* misrepresented as being the original software.
* 3. This notice may not be removed or altered from any source distribution.
*/

#ifndef B2_DYNAMIC_TREE_H
#define B2_DYNAMIC_TREE_H

#include "clam.h"
#include <cassert>
#include <cstdint>
#include <stack>
#include <vector>

#define b2_nullNode (-1)

struct AABB{
    double area(void)const{
        auto size = max_ - min_;
        return 2.0 * (size[0] * size[1] + size[0] * size[2] + size[1] * size[2]);
    }

    void combine(const AABB& a, const AABB& b){
        *this = a;
        for(int i = 0; i < 3; ++i){
            if(b.min_[i] < min_[i]) min_[i] = b.min_[i];
            if(b.max_[i] > max_[i]) max_[i] = b.max_[i];
        }
    }

    void translate(const clam::Vec3d& dr){
        min_ += dr;
        max_ += dr;
    }

    void fatten(double value){
        max_ += value;
        min_ -= value;
    }

    clam::Vec3d min_;
    clam::Vec3d max_;
};

inline bool aabb_overlap(const AABB& a, const AABB& b){
    for(size_t d = 0; d < 3; ++d){
        if(a.min_[d] > b.max_[d] || a.max_[d] < b.min_[d]) return false;
    }
    return true;
}

/// A node in the dynamic tree. The client does not interact with this directly.
struct TreeNode
{
	bool IsLeaf() const
	{
		return child1 == b2_nullNode;
	}

	/// Enlarged AABB
	AABB aabb;

	void* userData;

	union
	{
		int32_t parent;
		int32_t next;
	};

	int32_t child1;
	int32_t child2;

	// leaf = 0, free node = -1
	int32_t height;
};

/// A dynamic AABB tree broad-phase, inspired by Nathanael Presson's btDbvt.
/// A dynamic tree arranges data in a binary tree to accelerate
/// queries such as volume queries and ray casts. Leafs are proxies
/// with an AABB. In the tree we expand the proxy AABB by b2_fatAABBFactor
/// so that the proxy AABB is bigger than the client object. This allows the client
/// object to move by small amounts without triggering a tree update.
///
/// Nodes are pooled and relocatable, so we use node indices rather than pointers.
class AABBTree
{
public:
	/// Constructing the tree initializes the node pool.
	AABBTree();

	/// Destroy the tree, freeing the node pool.
	~AABBTree();

	/// Create a proxy. Provide a tight fitting AABB and a userData pointer.
	int32_t CreateProxy(const AABB& aabb, void* userData);

	/// Destroy a proxy. This asserts if the id is invalid.
	void DestroyProxy(int32_t proxyId);

	/// Move a proxy with a swepted AABB. If the proxy has moved outside of its fattened AABB,
	/// then the proxy is removed from the tree and re-inserted. Otherwise
	/// the function returns immediately.
	/// @return true if the proxy was re-inserted.
	bool MoveProxy(int32_t proxyId, const AABB& aabb1);

	/// Get proxy user data.
	/// @return the proxy user data or 0 if the id is invalid.
	void* GetUserData(int32_t proxyId) const;

	/// Get the fat AABB for a proxy.
	const AABB& GetFatAABB(int32_t proxyId) const;

	/// Query an AABB for overlapping proxies. The callback class
	/// is called for each proxy that overlaps the supplied AABB.
	template <typename T>
	void Query(const T& callback, const AABB& aabb) const;

	/// Validate this tree. For testing.
	void Validate() const;

	/// Compute the height of the binary tree in O(N) time. Should not be
	/// called often.
	int32_t GetHeight() const;

	/// Get the maximum balance of an node in the tree. The balance is the difference
	/// in height of the two children of a node.
	int32_t GetMaxBalance() const;

	/// Get the ratio of the sum of the node areas to the root area.
	double GetAreaRatio() const;

	/// Build an optimal tree. Very expensive. For testing.
	void RebuildBottomUp();

	/// Shift the world origin. Useful for large worlds.
	/// The shift formula is: position -= newOrigin
	/// @param newOrigin the new origin with respect to the old origin
	void ShiftOrigin(const clam::Vec3d& newOrigin);

private:

	int32_t AllocateNode();
	void FreeNode(int32_t node);

	void InsertLeaf(int32_t node);
	void RemoveLeaf(int32_t node);

	int32_t Balance(int32_t index);

	int32_t ComputeHeight() const;
	int32_t ComputeHeight(int32_t nodeId) const;

	void ValidateStructure(int32_t index) const;
	void ValidateMetrics(int32_t index) const;

	int32_t m_root;

	TreeNode* m_nodes;
	int32_t m_nodeCount;
	int32_t m_nodeCapacity;

	int32_t m_freeList;

	/// This is used to incrementally traverse the tree for re-balancing.
	uint32_t m_path;

	int32_t m_insertionCount;
};

inline void* AABBTree::GetUserData(int32_t proxyId) const
{
	assert(0 <= proxyId && proxyId < m_nodeCapacity);
	return m_nodes[proxyId].userData;
}

inline const AABB& AABBTree::GetFatAABB(int32_t proxyId) const
{
	assert(0 <= proxyId && proxyId < m_nodeCapacity);
	return m_nodes[proxyId].aabb;
}

template <typename T>
inline void AABBTree::Query(const T& callback, const AABB& aabb) const
{
    std::stack<int32_t, std::vector<int32_t>> stack;
	stack.push(m_root);

	while (stack.size() > 0)
	{
		int32_t nodeId = stack.top();
        stack.pop();
		if (nodeId == b2_nullNode)
		{
			continue;
		}

		const TreeNode* node = m_nodes + nodeId;

		if (aabb_overlap(node->aabb, aabb))
		{
			if (node->IsLeaf())
			{
				bool proceed = callback(nodeId);
				if (proceed == false)
				{
					return;
				}
			}
			else
			{
				stack.push(node->child1);
				stack.push(node->child2);
			}
		}
	}
}

#endif
