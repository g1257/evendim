#ifndef GRAPH_HH
#define GRAPH_HH
#include "PsimagLite.h"
#include "InputCheck.h"
#include "InputNg.h"
#include "BitManip.h"

namespace Gep {

class Graph {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	using LongUintType = long unsigned int;

	Graph(PsimagLite::String graphFile, SizeType vertices = 0, bool periodic = false)
	    : graphFile_(graphFile), vertices_(vertices), isConnected_(false)
	{
		if (graphFile_ == "zz") {
			createChain(periodic);
			return;
		}

		const bool fromFile = (graphFile_.substr(0, 5) == filePrefix());
		if (!fromFile)
			err("Unknown named graph " + graphFile + "\n");

		graphFile_ = graphFile_.substr(5, graphFile_.length() - filePrefix().length());

		PsimagLite::String data;
		PsimagLite::InputNg<InputCheck>::Writeable::readFile(data, graphFile_);

		data = discardComments(data, '#');
		stripLeadingChars(data);
		stripTrailingChars(data);
		loadGraphFromString(data);
		LongUintType state = getState();
		isConnected_ = isConnectedRunOnce(state);
	}

	Graph(LongUintType state, SizeType vertices)
	    : vertices_(vertices), isConnected_(false)
	{
		if (vertices < 2)
			err("Graph::ctor(): cannot construct Graph with less than two vertices\n");

		const SizeType pyramid = ((vertices_ - 1)*vertices_)/2;
		for (SizeType site1 = 0; site1 < vertices - 1; ++site1) {
			SizeType offset1 = findOffset(site1, pyramid);
			for (SizeType site2 = site1 + 1 ; site2 < vertices; ++site2) {
				const SizeType offset12 = offset1 + site2;
				const LongUintType mask = (1<<offset12);
				if ((state & mask) == 0) continue;
			}
		}

		isConnected_ = isConnectedRunOnce(state);
	}

	static PsimagLite::String filePrefix()
	{
		return "file:";
	}

	SizeType vertices() const
	{
		return vertices_;
	}

	const VectorSizeType& neighbors(SizeType site) const
	{
		assert(site < allNeighbors_.size());
		return allNeighbors_[site];
	}

	bool isConnected() const { return isConnected_; }

	friend std::ostream& operator<<(std::ostream& os, const Graph& graph)
	{
		const SizeType n = graph.vertices();
		if (n < 2) err("Cannot print a graph with less than two vertices\n");

		for (SizeType i = 0; i < n - 1; ++i) {
			PsimagLite::String str;
			graph.qaoaForVertex(str, i);
			os<<str<<"\n";
		}

		return os;
	}

private:

	void qaoaForVertex(PsimagLite::String& str, SizeType vertex) const
	{
		assert(vertex < allNeighbors_.size());
		neighborsToQaoa(str, vertex, allNeighbors_[vertex]);
	}

	void neighborsToQaoa(PsimagLite::String& str,
	                     SizeType vertex,
	                     const VectorSizeType& v) const
	{
		for (SizeType i = vertex + 1; i < vertices_; ++i) {
			const unsigned char c = (std::find(v.begin(), v.end(), i) == v.end()) ? '0' : '1';
			str += c;
		}
	}

	void neighborsToQaoa(LongUintType& state,
	                     SizeType& location,
	                     SizeType vertex,
	                     const VectorSizeType& v) const
	{

		for (SizeType i = vertex + 1; i < vertices_; ++i) {
			if (std::find(v.begin(), v.end(), i) == v.end()) continue;
			const LongUintType mask = (1<<location);
			state |= mask;
			checkLocation(location);
			++location;
		}
	}

	void createChain(bool periodic)
	{
		for (SizeType vertex = 0; vertex < vertices_; ++vertex) {
			SizeType nextSite = vertex + 1;
			if (nextSite == vertices_) {
				if (!periodic) {
					allNeighbors_.push_back(VectorSizeType());
					break;
				}

				nextSite = 0;
			}

			VectorSizeType tmpVector(1, nextSite);
			allNeighbors_.push_back(tmpVector);
		}
	}

	// Graph:1.0:Adjacency
	// #This is a comment
	// number of vertices here;
	// 0: 1;
	// 1: 2;
	// 2: 3;
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

		ind = readUntil(str, ind, data, ';');

		SizeType vertices = PsimagLite::atoi(str);
		if (vertices_ > 0 && vertices_ != vertices)
			err("Expected " + ttos(vertices_) + " in " + graphFile_ +
			    ", but got " + str + " instead.\n");

		vertices_ = vertices;

		while (!isTheEnd(ind, data)) {
			ind = readUntil(str, ind, data, ':');
			str = stripBlanks(str);
			SizeType site = PsimagLite::atoi(str);

			ind = readUntil(str, ind, data, ';');
			tokens.clear();
			PsimagLite::split(tokens, str, " ");
			if (tokens.size() == 0)
				err("No connections given for listed site " + ttos(site) + "\n");

			if (tokens.size() + 1 >= vertices_)
				err("Too many connections given for listed site " + ttos(site) + "\n");

			addToNeighbors(site, tokens);
		}
	}

	void addToNeighbors(SizeType site, const VectorStringType& neighs)
	{
		VectorSizeType tmpVector(neighs.size());
		for (auto it = neighs.begin(); it != neighs.end(); ++it) {
			const SizeType site2 = PsimagLite::atoi(*it);
			if (std::find(tmpVector.begin(), tmpVector.end(), site2) != tmpVector.end())
				err("Site already added to neighbors of " + ttos(site) + "\n");
			tmpVector[it - neighs.begin()] = PsimagLite::atoi(*it);
		}

		if (allNeighbors_.size() == 0)
			allNeighbors_.resize(vertices_);

		assert(site < allNeighbors_.size());
		allNeighbors_[site] = tmpVector;
	}

	void loadFromGraphQaoa(PsimagLite::String data, SizeType ind, PsimagLite::String str)
	{
		vertices_ = readOrderGraphQaoa(str);
		if (vertices_ < 2) err("loadFromGraphQaoa: Only one vertex found!?\n");

		for (SizeType i = 0; i < vertices_ - 1; ++i) {
			ind = readUntil(str, ind, data, '\n');
			procVertexGraphQaoa(i, str);
		}
	}

	void procVertexGraphQaoa(SizeType site, PsimagLite::String str)
	{
		const SizeType total = str.size();
		VectorSizeType tmpVector;
		for (SizeType pos = 0; pos < total; ++pos) {
			const unsigned char c = str[pos];
			if (c == '0') continue;
			if (c != '1') err("procVertexGraphQaoa: Expected 0 or 1\n");
			tmpVector.push_back(pos + site + 1);
		}

		if (allNeighbors_.size() == 0)
			allNeighbors_.resize(vertices_);

		assert(site < allNeighbors_.size());
		allNeighbors_[site] = tmpVector;
	}

	SizeType findOffset(SizeType site1, SizeType pyramid)
	{
		const SizeType tmp1 = vertices_ - site1 - 1;
		const SizeType tmp2 = tmp1*(tmp1 + 1);
		const SizeType tmp3 = tmp2/2;
		assert(pyramid >= tmp3);
		return pyramid - tmp3;
	}

	bool isConnectedRunOnce(LongUintType state) const
	{
		if (PsimagLite::BitManip::countKernighan(state) + 1 < vertices_) return false;
		return true;
	}

	LongUintType getState() const
	{
		LongUintType state = 0;
		SizeType location = 0;
		for (SizeType vertex = 0; vertex < vertices_ - 1; ++vertex) {
			neighborsToQaoa(state, location, vertex, allNeighbors_[vertex]);
		}

		return state;
	}

	void checkLocation(SizeType location) const
	{
		assert(location < ((vertices_ - 1)*vertices_)/2);
	}

	static SizeType readOrderGraphQaoa(PsimagLite::String str)
	{
		const SizeType total = str.size();
		PsimagLite::String buffer;
		SizeType ind = total;
		for (SizeType pos = 0; pos < total; ++pos) {
			ind = total - pos - 1;
			if (str[ind] == '.') continue;
			if (!std::isdigit(str[ind])) break;
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
			if (c2 == c) break;
			if (c2 == '\n') c2 = ' ';
			buffer += c2;
		}

		while (data[++pos] == c) {}
		return pos;
	}

	static PsimagLite::String discardComments(PsimagLite::String data, unsigned char c)
	{
		PsimagLite::String newData;
		SizeType ind = 0;
		while (ind < data.size()) {
			PsimagLite::String buffer;
			ind = readUntil(buffer, ind, data, '\n');
			if (buffer.size() > 0 && buffer[0] == c) continue;
			newData += buffer + "\n";
		}

		return newData;
	}

	static void stripLeadingChars(PsimagLite::String& data)
	{
		const SizeType total = data.size();
		SizeType pos = 0;
		for (; pos < total; ++pos) {
			if (!isBlankChar(data[pos])) break;
		}

		data = data.substr(pos, total - pos);
	}

	static void stripTrailingChars(PsimagLite::String& data)
	{
		const SizeType total = data.size();
		SizeType ind = total;
		for (SizeType pos = 0; pos < total; ++pos) {
			ind = total - pos - 1;
			if (!isBlankChar(data[ind])) break;
		}

		data = data.substr(0, ind + 2);
		data[ind + 1] = '\n';
	}

	static bool isTheEnd(SizeType ind, PsimagLite::String data)
	{
		if (ind >= data.size()) return true;
		return (ind + 1 == data.size() && data[ind] == '\n');
	}

	static PsimagLite::String stripBlanks(PsimagLite::String str)
	{
		const SizeType total = str.size();
		PsimagLite::String buffer;
		for (SizeType pos = 0; pos < total; ++pos) {
			if (isBlankChar(str[pos])) continue;
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
	VectorVectorSizeType allNeighbors_;
};
}
#endif // GRAPH_HH
