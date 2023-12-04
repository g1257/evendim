#ifndef GRAPH_HH
#define GRAPH_HH
#include "../Engine/InputCheck.h"
#include "BitManip.h"
#include "InputNg.h"
#include "PsimagLite.h"

namespace Gep {

class Graph {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	using VectorBoolType = PsimagLite::Vector<bool>::Type;
	using VectorVectorBoolType = PsimagLite::Vector<VectorBoolType>::Type;
	typedef PsimagLite::Vector<PsimagLite::String>::Type
	    VectorStringType;
	using LongUintType = long unsigned int;

	Graph(PsimagLite::String graphFile, SizeType vertices = 0, bool periodic = false)
	    : graphFile_(graphFile)
	    , vertices_(vertices)
	    , isConnected_(false)
	{
		if (graphFile_ == "zz") {
			createChain(periodic);
			return;
		}

		const bool fromFile = (graphFile_.substr(0, 5) == filePrefix());
		if (!fromFile)
			err("Unknown named graph " + graphFile + "\n");

		graphFile_ = graphFile_.substr(
		    5, graphFile_.length() - filePrefix().length());

		PsimagLite::String data;
		PsimagLite::InputNg<InputCheck>::Writeable::readFile(
		    data, graphFile_);

		data = discardComments(data, '#');
		stripLeadingChars(data);
		stripTrailingChars(data);
		loadGraphFromString(data);
		LongUintType state = getState();
		isConnected_ = isConnectedRunOnce(state);
	}

	Graph(LongUintType state, SizeType vertices)
	    : vertices_(vertices)
	    , isConnected_(false)
	{
		if (vertices < 2)
			err("Graph::ctor(): cannot construct Graph "
			    "with less "
			    "than two vertices\n");

		triangular_.resize(vertices_ - 1);
		const SizeType pyramid = ((vertices_ - 1) * vertices_) / 2;
		for (SizeType site1 = 0; site1 < vertices - 1;
		     ++site1) {
			SizeType offset1 = findOffset(site1, pyramid);
			VectorBoolType tmpVector(vertices - site1 - 1,
			                         false);
			for (SizeType site2 = site1 + 1;
			     site2 < vertices;
			     ++site2) {
				const SizeType j = site2 - site1 - 1;
				const SizeType offset12 = offset1 + j;
				const LongUintType mask = (1 << offset12);
				tmpVector[j] = (state & mask) ? 1 : 0;
			}

			triangular_[site1] = tmpVector;
		}

		isConnected_ = isConnectedRunOnce(state);
	}

	static PsimagLite::String filePrefix() { return "file:"; }

	SizeType vertices() const { return vertices_; }

	bool isConnected() const { return isConnected_; }

	bool connected(SizeType site1, SizeType site2) const
	{
		assert(site1 < vertices_ && site2 < vertices_);
		if (site1 == site2)
			return false;

		const SizeType minSite = (site1 < site2) ? site1 : site2;
		const SizeType maxSite = (site1 < site2) ? site2 : site1;
		assert(minSite < triangular_.size());
		assert(minSite < maxSite);
		const SizeType diff = maxSite - minSite;
		assert(diff > 0);
		assert(triangular_[minSite].size() + 1 > diff);
		return triangular_[minSite][diff - 1];
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const Graph& graph)
	{
		const SizeType n = graph.vertices();
		if (n < 2)
			err("Cannot print a graph with less than two "
			    "vertices\n");

		for (SizeType i = 0; i < n - 1; ++i) {
			PsimagLite::String str;
			graph.qaoaForVertex(str, i);
			os << str << "\n";
		}

		return os;
	}

private:

	void qaoaForVertex(PsimagLite::String& str,
	                   SizeType vertex) const
	{
		assert(vertex < triangular_.size());
		neighborsToQaoa(str, triangular_[vertex]);
	}

	void neighborsToQaoa(PsimagLite::String& str,
	                     const VectorBoolType& v) const
	{
		for (SizeType i = 0; i < v.size(); ++i) {
			const unsigned char c = (v[i]) ? '1' : '0';
			str += c;
		}
	}

	void neighborsToQaoa(LongUintType& state, SizeType& location, const VectorBoolType& v) const
	{
		const SizeType n = v.size();
		for (SizeType i = 0; i < n; ++i) {
			if (!v[i]) {
				++location;
				continue;
			}

			const LongUintType mask = (1 << location);
			state |= mask;
			checkLocation(location);
			++location;
		}
	}

	void createChain(bool periodic)
	{
		assert(vertices_ > 1);
		for (SizeType vertex = 0; vertex < vertices_ - 1;
		     ++vertex) {
			VectorBoolType v(vertices_ - vertex - 1, false);
			v[0] = true;
			if (periodic && vertex == 0 && v.size() >= 2) {
				assert(v.size() > 1);
				v[v.size() - 1] = true;
			}

			triangular_.push_back(v);
		}
	}

	// Graph:1.0:Adjacency
	// #This is a comment
	void loadGraphFromString(PsimagLite::String data)
	{
		static const PsimagLite::String GRAPH = "Graph";

		PsimagLite::String str;
		SizeType ind = readUntil(str, 0, data, '\n');

		VectorStringType tokens;
		PsimagLite::split(tokens, str, ":");
		if (tokens.size() == 1 && tokens[0].substr(0, GRAPH.length()) == GRAPH) {
			loadFromGraphQaoa(data, ind, str);
			return;
		}

		if (tokens.size() < 2 || tokens[0] != "Graph")
			err("Expected Graph:VersionNumber not " + str + " in " + graphFile_ + "\n");

		ind = readUntil(str, ind, data, '\n');

		SizeType vertices = PsimagLite::atoi(str);
		if (vertices_ > 0 && vertices_ != vertices)
			err("Expected " + ttos(vertices_) + " in " + graphFile_ + ", but got " + str + " instead.\n");

		vertices_ = vertices;

		SizeType count = vertices_ - 1;
		SizeType site = 0;
		while (!isTheEnd(ind, data)) {
			ind = readUntil(str, ind, data, '\n');

			assert(count > 0);
			if (str.size() != count)
				err("Expected " + ttos(count) + " numbers for site " + ttos(site) + ", not " + ttos(str.size()) + "\n");

			addToNeighbors(site, str);
			--count;
			++site;
		}
	}

	void addToNeighbors(const SizeType site,
	                    PsimagLite::String neighs)
	{
		assert(site + 1 < vertices_);
		assert(neighs.size() == vertices_ - site - 1);
		VectorBoolType tmpVector(neighs.size());
		for (SizeType i = site + 1; i < vertices_; ++i) {
			const SizeType j = i - site - 1;
			const unsigned char c = neighs[j];
			if (c != '0' && c != '1')
				err("addToNeighbors: adjancency matrix "
				    "found "
				    + neighs + " not 0 or 1\n");

			assert(j < tmpVector.size());
			tmpVector[j] = (c == '0') ? false : true;
		}

		if (triangular_.size() == 0)
			triangular_.resize(vertices_ - 1);

		assert(site < triangular_.size());
		triangular_[site] = tmpVector;
	}

	void loadFromGraphQaoa(PsimagLite::String data, SizeType ind, PsimagLite::String str)
	{
		vertices_ = readOrderGraphQaoa(str);
		if (vertices_ < 2)
			err("loadFromGraphQaoa: Only one vertex "
			    "found!?\n");

		for (SizeType i = 0; i < vertices_ - 1; ++i) {
			ind = readUntil(str, ind, data, '\n');
			addToNeighbors(i, str);
		}
	}

	SizeType findOffset(SizeType site1, SizeType pyramid)
	{
		const SizeType tmp1 = vertices_ - site1 - 1;
		const SizeType tmp2 = tmp1 * (tmp1 + 1);
		const SizeType tmp3 = tmp2 / 2;
		assert(pyramid >= tmp3);
		return pyramid - tmp3;
	}

	bool isConnectedRunOnce(LongUintType state) const
	{
		if (PsimagLite::BitManip::countKernighan(state) + 1 < vertices_)
			return false;
		VectorBoolType visited(vertices_, false);
		visitVertex(visited, 0);

		return allAreTrue(visited);
	}

	void visitVertex(VectorBoolType& visited, SizeType vertex) const
	{
		assert(vertex < visited.size());
		if (visited[vertex])
			return;
		visited[vertex] = true;
		if (vertex + 1 == vertices_)
			return;
		assert(vertex < triangular_.size());
		for (SizeType i = 0; i < triangular_[vertex].size();
		     ++i) {
			if (!triangular_[vertex][i])
				continue;
			const SizeType vertex2 = vertex + i + 1;
			visitVertex(visited, vertex2);
		}
	}

	static bool allAreTrue(const VectorBoolType& v)
	{
		for (SizeType i = 0; i < v.size(); ++i)
			if (!v[i])
				return false;

		return true;
	}

	LongUintType getState() const
	{
		LongUintType state = 0;
		SizeType location = 0;
		for (SizeType vertex = 0; vertex < vertices_ - 1;
		     ++vertex) {
			neighborsToQaoa(state, location, triangular_[vertex]);
		}

		return state;
	}

	void checkLocation(SizeType location) const
	{
		assert(location < ((vertices_ - 1) * vertices_) / 2);
	}

	static SizeType readOrderGraphQaoa(PsimagLite::String str)
	{
		const SizeType total = str.size();
		PsimagLite::String buffer;
		SizeType ind = total;
		for (SizeType pos = 0; pos < total; ++pos) {
			ind = total - pos - 1;
			if (str[ind] == '.')
				continue;
			if (!std::isdigit(str[ind]))
				break;
			buffer += str[ind];
		}

		const SizeType total2 = buffer.size();
		str = "";
		for (SizeType pos = 0; pos < total2; ++pos) {
			const SizeType ind2 = total2 - pos - 1;
			str += buffer[ind2];
		}

		return PsimagLite::atoi(str);
	}

	static SizeType readUntil(PsimagLite::String& buffer,
	                          SizeType ind,
	                          PsimagLite::String data,
	                          unsigned char c)
	{
		buffer = "";

		const SizeType total = data.size();
		SizeType pos = ind;
		for (; pos < total; ++pos) {
			unsigned char c2 = data[pos];
			if (c2 == c)
				break;
			if (c2 == '\n')
				c2 = ' ';
			buffer += c2;
		}

		while (data[++pos] == c) {
		}
		return pos;
	}

	static PsimagLite::String
	discardComments(PsimagLite::String data, unsigned char c)
	{
		PsimagLite::String newData;
		SizeType ind = 0;
		while (ind < data.size()) {
			PsimagLite::String buffer;
			ind = readUntil(buffer, ind, data, '\n');
			if (buffer.size() > 0 && buffer[0] == c)
				continue;
			newData += buffer + "\n";
		}

		return newData;
	}

	static void stripLeadingChars(PsimagLite::String& data)
	{
		const SizeType total = data.size();
		SizeType pos = 0;
		for (; pos < total; ++pos) {
			if (!isBlankChar(data[pos]))
				break;
		}

		data = data.substr(pos, total - pos);
	}

	static void stripTrailingChars(PsimagLite::String& data)
	{
		const SizeType total = data.size();
		SizeType ind = total;
		for (SizeType pos = 0; pos < total; ++pos) {
			ind = total - pos - 1;
			if (!isBlankChar(data[ind]))
				break;
		}

		data = data.substr(0, ind + 2);
		data[ind + 1] = '\n';
	}

	static bool isTheEnd(SizeType ind, PsimagLite::String data)
	{
		if (ind >= data.size())
			return true;
		return (ind + 1 == data.size() && data[ind] == '\n');
	}

	static PsimagLite::String stripBlanks(PsimagLite::String str)
	{
		const SizeType total = str.size();
		PsimagLite::String buffer;
		for (SizeType pos = 0; pos < total; ++pos) {
			if (isBlankChar(str[pos]))
				continue;
			buffer += str[pos];
		}

		return buffer;
	}

	static bool isBlankChar(unsigned char c)
	{
		return (c == ' ' || c == '\n' || c == '\t');
	}

	PsimagLite::String graphFile_;
	SizeType vertices_;
	bool isConnected_;
	VectorVectorBoolType triangular_;
};
} // namespace Gep
#endif // GRAPH_HH
