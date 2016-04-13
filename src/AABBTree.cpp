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

#include "AABBTree.h"
#include <cstring>
#include <cstdlib>
#include <algorithm>

AABBTree::AABBTree()
{
	m_root = b2_nullNode;

	m_nodeCapacity = 16;
	m_nodeCount = 0;
	m_nodes = (TreeNode*)malloc(m_nodeCapacity * sizeof(TreeNode));
	memset(m_nodes, 0, m_nodeCapacity * sizeof(TreeNode));

	// Build a linked list for the free list.
	for (int32_t i = 0; i < m_nodeCapacity - 1; ++i)
	{
		m_nodes[i].next = i + 1;
		m_nodes[i].height = -1;
	}
	m_nodes[m_nodeCapacity-1].next = b2_nullNode;
	m_nodes[m_nodeCapacity-1].height = -1;
	m_freeList = 0;

	m_path = 0;

	m_insertionCount = 0;
}

AABBTree::~AABBTree()
{
	// This frees the entire tree in one shot.
	free(m_nodes);
}

// Allocate a node from the pool. Grow the pool if necessary.
int32_t AABBTree::AllocateNode()
{
	// Expand the node pool as needed.
	if (m_freeList == b2_nullNode)
	{
		assert(m_nodeCount == m_nodeCapacity);

		// The free list is empty. Rebuild a bigger pool.
		TreeNode* oldNodes = m_nodes;
		m_nodeCapacity *= 2;
		m_nodes = (TreeNode*)malloc(m_nodeCapacity * sizeof(TreeNode));
		memcpy(m_nodes, oldNodes, m_nodeCount * sizeof(TreeNode));
		free(oldNodes);

		// Build a linked list for the free list. The parent
		// pointer becomes the "next" pointer.
		for (int32_t i = m_nodeCount; i < m_nodeCapacity - 1; ++i)
		{
			m_nodes[i].next = i + 1;
			m_nodes[i].height = -1;
		}
		m_nodes[m_nodeCapacity-1].next = b2_nullNode;
		m_nodes[m_nodeCapacity-1].height = -1;
		m_freeList = m_nodeCount;
	}

	// Peel a node off the free list.
	int32_t nodeId = m_freeList;
	m_freeList = m_nodes[nodeId].next;
	m_nodes[nodeId].parent = b2_nullNode;
	m_nodes[nodeId].child1 = b2_nullNode;
	m_nodes[nodeId].child2 = b2_nullNode;
	m_nodes[nodeId].height = 0;
	m_nodes[nodeId].userData = NULL;
	++m_nodeCount;
	return nodeId;
}

// Return a node to the pool.
void AABBTree::FreeNode(int32_t nodeId)
{
	assert(0 <= nodeId && nodeId < m_nodeCapacity);
	assert(0 < m_nodeCount);
	m_nodes[nodeId].next = m_freeList;
	m_nodes[nodeId].height = -1;
	m_freeList = nodeId;
	--m_nodeCount;
}

// Create a proxy in the tree as a leaf node. We return the index
// of the node instead of a pointer so that we can grow
// the node pool.
int32_t AABBTree::CreateProxy(const AABB& aabb, void* userData)
{
	int32_t proxyId = AllocateNode();

	// Fatten the aabb.
	m_nodes[proxyId].aabb.min_ = aabb.min_;
	m_nodes[proxyId].aabb.max_ = aabb.max_;
	m_nodes[proxyId].userData = userData;
	m_nodes[proxyId].height = 0;

	InsertLeaf(proxyId);

	return proxyId;
}

void AABBTree::DestroyProxy(int32_t proxyId)
{
	assert(0 <= proxyId && proxyId < m_nodeCapacity);
	assert(m_nodes[proxyId].IsLeaf());

	RemoveLeaf(proxyId);
	FreeNode(proxyId);
}

bool AABBTree::MoveProxy(int32_t proxyId, const AABB& aabb)
{
	assert(0 <= proxyId && proxyId < m_nodeCapacity);

	assert(m_nodes[proxyId].IsLeaf());

	RemoveLeaf(proxyId);

	m_nodes[proxyId].aabb = aabb;

	InsertLeaf(proxyId);
	return true;
}

void AABBTree::InsertLeaf(int32_t leaf)
{
	++m_insertionCount;

	if (m_root == b2_nullNode)
	{
		m_root = leaf;
		m_nodes[m_root].parent = b2_nullNode;
		return;
	}

	// Find the best sibling for this node
	AABB leafAABB = m_nodes[leaf].aabb;
	int32_t index = m_root;
	while (m_nodes[index].IsLeaf() == false)
	{
		int32_t child1 = m_nodes[index].child1;
		int32_t child2 = m_nodes[index].child2;

		double area = m_nodes[index].aabb.area();

		AABB combinedAABB;
		combinedAABB.combine(m_nodes[index].aabb, leafAABB);
		double combinedArea = combinedAABB.area();

		// Cost of creating a new parent for this node and the new leaf
		double cost = 2.0f * combinedArea;

		// Minimum cost of pushing the leaf further down the tree
		double inheritanceCost = 2.0f * (combinedArea - area);

		// Cost of descending into child1
		double cost1;
		if (m_nodes[child1].IsLeaf())
		{
			AABB aabb;
			aabb.combine(leafAABB, m_nodes[child1].aabb);
			cost1 = aabb.area() + inheritanceCost;
		}
		else
		{
			AABB aabb;
			aabb.combine(leafAABB, m_nodes[child1].aabb);
			double oldArea = m_nodes[child1].aabb.area();
			double newArea = aabb.area();
			cost1 = (newArea - oldArea) + inheritanceCost;
		}

		// Cost of descending into child2
		double cost2;
		if (m_nodes[child2].IsLeaf())
		{
			AABB aabb;
			aabb.combine(leafAABB, m_nodes[child2].aabb);
			cost2 = aabb.area() + inheritanceCost;
		}
		else
		{
			AABB aabb;
			aabb.combine(leafAABB, m_nodes[child2].aabb);
			double oldArea = m_nodes[child2].aabb.area();
			double newArea = aabb.area();
			cost2 = newArea - oldArea + inheritanceCost;
		}

		// Descend according to the minimum cost.
		if (cost < cost1 && cost < cost2)
		{
			break;
		}

		// Descend
		if (cost1 < cost2)
		{
			index = child1;
		}
		else
		{
			index = child2;
		}
	}

	int32_t sibling = index;

	// Create a new parent.
	int32_t oldParent = m_nodes[sibling].parent;
	int32_t newParent = AllocateNode();
	m_nodes[newParent].parent = oldParent;
	m_nodes[newParent].userData = NULL;
	m_nodes[newParent].aabb.combine(leafAABB, m_nodes[sibling].aabb);
	m_nodes[newParent].height = m_nodes[sibling].height + 1;

	if (oldParent != b2_nullNode)
	{
		// The sibling was not the root.
		if (m_nodes[oldParent].child1 == sibling)
		{
			m_nodes[oldParent].child1 = newParent;
		}
		else
		{
			m_nodes[oldParent].child2 = newParent;
		}

		m_nodes[newParent].child1 = sibling;
		m_nodes[newParent].child2 = leaf;
		m_nodes[sibling].parent = newParent;
		m_nodes[leaf].parent = newParent;
	}
	else
	{
		// The sibling was the root.
		m_nodes[newParent].child1 = sibling;
		m_nodes[newParent].child2 = leaf;
		m_nodes[sibling].parent = newParent;
		m_nodes[leaf].parent = newParent;
		m_root = newParent;
	}

	// Walk back up the tree fixing heights and AABBs
	index = m_nodes[leaf].parent;
	while (index != b2_nullNode)
	{
		index = Balance(index);

		int32_t child1 = m_nodes[index].child1;
		int32_t child2 = m_nodes[index].child2;

		assert(child1 != b2_nullNode);
		assert(child2 != b2_nullNode);

		m_nodes[index].height = 1 + std::max(m_nodes[child1].height, m_nodes[child2].height);
		m_nodes[index].aabb.combine(m_nodes[child1].aabb, m_nodes[child2].aabb);

		index = m_nodes[index].parent;
	}

	//Validate();
}

void AABBTree::RemoveLeaf(int32_t leaf)
{
	if (leaf == m_root)
	{
		m_root = b2_nullNode;
		return;
	}

	int32_t parent = m_nodes[leaf].parent;
	int32_t grandParent = m_nodes[parent].parent;
	int32_t sibling;
	if (m_nodes[parent].child1 == leaf)
	{
		sibling = m_nodes[parent].child2;
	}
	else
	{
		sibling = m_nodes[parent].child1;
	}

	if (grandParent != b2_nullNode)
	{
		// Destroy parent and connect sibling to grandParent.
		if (m_nodes[grandParent].child1 == parent)
		{
			m_nodes[grandParent].child1 = sibling;
		}
		else
		{
			m_nodes[grandParent].child2 = sibling;
		}
		m_nodes[sibling].parent = grandParent;
		FreeNode(parent);

		// Adjust ancestor bounds.
		int32_t index = grandParent;
		while (index != b2_nullNode)
		{
			index = Balance(index);

			int32_t child1 = m_nodes[index].child1;
			int32_t child2 = m_nodes[index].child2;

			m_nodes[index].aabb.combine(m_nodes[child1].aabb, m_nodes[child2].aabb);
			m_nodes[index].height = 1 + std::max(m_nodes[child1].height, m_nodes[child2].height);

			index = m_nodes[index].parent;
		}
	}
	else
	{
		m_root = sibling;
		m_nodes[sibling].parent = b2_nullNode;
		FreeNode(parent);
	}

	//Validate();
}

// Perform a left or right rotation if node A is imbalanced.
// Returns the new root index.
int32_t AABBTree::Balance(int32_t iA)
{
	assert(iA != b2_nullNode);

	TreeNode* A = m_nodes + iA;
	if (A->IsLeaf() || A->height < 2)
	{
		return iA;
	}

	int32_t iB = A->child1;
	int32_t iC = A->child2;
	assert(0 <= iB && iB < m_nodeCapacity);
	assert(0 <= iC && iC < m_nodeCapacity);

	TreeNode* B = m_nodes + iB;
	TreeNode* C = m_nodes + iC;

	int32_t balance = C->height - B->height;

	// Rotate C up
	if (balance > 1)
	{
		int32_t iF = C->child1;
		int32_t iG = C->child2;
		TreeNode* F = m_nodes + iF;
		TreeNode* G = m_nodes + iG;
		assert(0 <= iF && iF < m_nodeCapacity);
		assert(0 <= iG && iG < m_nodeCapacity);

		// Swap A and C
		C->child1 = iA;
		C->parent = A->parent;
		A->parent = iC;

		// A's old parent should point to C
		if (C->parent != b2_nullNode)
		{
			if (m_nodes[C->parent].child1 == iA)
			{
				m_nodes[C->parent].child1 = iC;
			}
			else
			{
				assert(m_nodes[C->parent].child2 == iA);
				m_nodes[C->parent].child2 = iC;
			}
		}
		else
		{
			m_root = iC;
		}

		// Rotate
		if (F->height > G->height)
		{
			C->child2 = iF;
			A->child2 = iG;
			G->parent = iA;
			A->aabb.combine(B->aabb, G->aabb);
			C->aabb.combine(A->aabb, F->aabb);

			A->height = 1 + std::max(B->height, G->height);
			C->height = 1 + std::max(A->height, F->height);
		}
		else
		{
			C->child2 = iG;
			A->child2 = iF;
			F->parent = iA;
			A->aabb.combine(B->aabb, F->aabb);
			C->aabb.combine(A->aabb, G->aabb);

			A->height = 1 + std::max(B->height, F->height);
			C->height = 1 + std::max(A->height, G->height);
		}

		return iC;
	}
	
	// Rotate B up
	if (balance < -1)
	{
		int32_t iD = B->child1;
		int32_t iE = B->child2;
		TreeNode* D = m_nodes + iD;
		TreeNode* E = m_nodes + iE;
		assert(0 <= iD && iD < m_nodeCapacity);
		assert(0 <= iE && iE < m_nodeCapacity);

		// Swap A and B
		B->child1 = iA;
		B->parent = A->parent;
		A->parent = iB;

		// A's old parent should point to B
		if (B->parent != b2_nullNode)
		{
			if (m_nodes[B->parent].child1 == iA)
			{
				m_nodes[B->parent].child1 = iB;
			}
			else
			{
				assert(m_nodes[B->parent].child2 == iA);
				m_nodes[B->parent].child2 = iB;
			}
		}
		else
		{
			m_root = iB;
		}

		// Rotate
		if (D->height > E->height)
		{
			B->child2 = iD;
			A->child1 = iE;
			E->parent = iA;
			A->aabb.combine(C->aabb, E->aabb);
			B->aabb.combine(A->aabb, D->aabb);

			A->height = 1 + std::max(C->height, E->height);
			B->height = 1 + std::max(A->height, D->height);
		}
		else
		{
			B->child2 = iE;
			A->child1 = iD;
			D->parent = iA;
			A->aabb.combine(C->aabb, D->aabb);
			B->aabb.combine(A->aabb, E->aabb);

			A->height = 1 + std::max(C->height, D->height);
			B->height = 1 + std::max(A->height, E->height);
		}

		return iB;
	}

	return iA;
}

int32_t AABBTree::GetHeight() const
{
	if (m_root == b2_nullNode)
	{
		return 0;
	}

	return m_nodes[m_root].height;
}

//
double AABBTree::GetAreaRatio() const
{
	if (m_root == b2_nullNode)
	{
		return 0.0f;
	}

	const TreeNode* root = m_nodes + m_root;
	double rootArea = root->aabb.area();

	double totalArea = 0.0f;
	for (int32_t i = 0; i < m_nodeCapacity; ++i)
	{
		const TreeNode* node = m_nodes + i;
		if (node->height < 0)
		{
			// Free node in pool
			continue;
		}

		totalArea += node->aabb.area();
	}

	return totalArea / rootArea;
}

// Compute the height of a sub-tree.
int32_t AABBTree::ComputeHeight(int32_t nodeId) const
{
	assert(0 <= nodeId && nodeId < m_nodeCapacity);
	TreeNode* node = m_nodes + nodeId;

	if (node->IsLeaf())
	{
		return 0;
	}

	int32_t height1 = ComputeHeight(node->child1);
	int32_t height2 = ComputeHeight(node->child2);
	return 1 + std::max(height1, height2);
}

int32_t AABBTree::ComputeHeight() const
{
	int32_t height = ComputeHeight(m_root);
	return height;
}

void AABBTree::ValidateStructure(int32_t index) const
{
	if (index == b2_nullNode)
	{
		return;
	}

	if (index == m_root)
	{
		assert(m_nodes[index].parent == b2_nullNode);
	}

	const TreeNode* node = m_nodes + index;

	int32_t child1 = node->child1;
	int32_t child2 = node->child2;

	if (node->IsLeaf())
	{
		assert(child1 == b2_nullNode);
		assert(child2 == b2_nullNode);
		assert(node->height == 0);
		return;
	}

	assert(0 <= child1 && child1 < m_nodeCapacity);
	assert(0 <= child2 && child2 < m_nodeCapacity);

	assert(m_nodes[child1].parent == index);
	assert(m_nodes[child2].parent == index);

	ValidateStructure(child1);
	ValidateStructure(child2);
}

void AABBTree::ValidateMetrics(int32_t index) const
{
	if (index == b2_nullNode)
	{
		return;
	}

	const TreeNode* node = m_nodes + index;

	int32_t child1 = node->child1;
	int32_t child2 = node->child2;

	if (node->IsLeaf())
	{
		assert(child1 == b2_nullNode);
		assert(child2 == b2_nullNode);
		assert(node->height == 0);
		return;
	}

	assert(0 <= child1 && child1 < m_nodeCapacity);
	assert(0 <= child2 && child2 < m_nodeCapacity);

	assert(node->height == 1 + std::max(m_nodes[child1].height, m_nodes[child2].height));

	AABB aabb;
	aabb.combine(m_nodes[child1].aabb, m_nodes[child2].aabb);

	assert(aabb.min_ == node->aabb.min_);
	assert(aabb.max_ == node->aabb.max_);

	ValidateMetrics(child1);
	ValidateMetrics(child2);
}

void AABBTree::Validate() const
{
	ValidateStructure(m_root);
	ValidateMetrics(m_root);

	int32_t freeCount = 0;
	int32_t freeIndex = m_freeList;
	while (freeIndex != b2_nullNode)
	{
		assert(0 <= freeIndex && freeIndex < m_nodeCapacity);
		freeIndex = m_nodes[freeIndex].next;
		++freeCount;
	}

	assert(GetHeight() == ComputeHeight());

	assert(m_nodeCount + freeCount == m_nodeCapacity);
}

int32_t AABBTree::GetMaxBalance() const
{
	int32_t maxBalance = 0;
	for (int32_t i = 0; i < m_nodeCapacity; ++i)
	{
		const TreeNode* node = m_nodes + i;
		if (node->height <= 1)
		{
			continue;
		}

		assert(node->IsLeaf() == false);

		int32_t child1 = node->child1;
		int32_t child2 = node->child2;
		int32_t balance = std::abs(m_nodes[child2].height - m_nodes[child1].height);
		maxBalance = std::max(maxBalance, balance);
	}

	return maxBalance;
}

void AABBTree::RebuildBottomUp()
{
	int32_t* nodes = (int32_t*)malloc(m_nodeCount * sizeof(int32_t));
	int32_t count = 0;

	// Build array of leaves. Free the rest.
	for (int32_t i = 0; i < m_nodeCapacity; ++i)
	{
		if (m_nodes[i].height < 0)
		{
			// free node in pool
			continue;
		}

		if (m_nodes[i].IsLeaf())
		{
			m_nodes[i].parent = b2_nullNode;
			nodes[count] = i;
			++count;
		}
		else
		{
			FreeNode(i);
		}
	}

	while (count > 1)
	{
		double minCost = std::numeric_limits<double>::max();
		int32_t iMin = -1, jMin = -1;
		for (int32_t i = 0; i < count; ++i)
		{
			AABB aabbi = m_nodes[nodes[i]].aabb;

			for (int32_t j = i + 1; j < count; ++j)
			{
				AABB aabbj = m_nodes[nodes[j]].aabb;
				AABB b;
				b.combine(aabbi, aabbj);
				double cost = b.area();
				if (cost < minCost)
				{
					iMin = i;
					jMin = j;
					minCost = cost;
				}
			}
		}

		int32_t index1 = nodes[iMin];
		int32_t index2 = nodes[jMin];
		TreeNode* child1 = m_nodes + index1;
		TreeNode* child2 = m_nodes + index2;

		int32_t parentIndex = AllocateNode();
		TreeNode* parent = m_nodes + parentIndex;
		parent->child1 = index1;
		parent->child2 = index2;
		parent->height = 1 + std::max(child1->height, child2->height);
		parent->aabb.combine(child1->aabb, child2->aabb);
		parent->parent = b2_nullNode;

		child1->parent = parentIndex;
		child2->parent = parentIndex;

		nodes[jMin] = nodes[count-1];
		nodes[iMin] = parentIndex;
		--count;
	}

	m_root = nodes[0];
	free(nodes);

	Validate();
}

void AABBTree::ShiftOrigin(const clam::Vec3d& newOrigin)
{
	// Build array of leaves. Free the rest.
	for (int32_t i = 0; i < m_nodeCapacity; ++i)
	{
		m_nodes[i].aabb.min_ -= newOrigin;
		m_nodes[i].aabb.max_ -= newOrigin;
	}
}
