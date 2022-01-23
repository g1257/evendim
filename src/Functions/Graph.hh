#ifndef GRAPH_HH
#define GRAPH_HH
#include "PsimagLite.h"
#include "InputCheck.h"
#include "InputNg.h"

namespace Gep {

class Graph {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	Graph(PsimagLite::String graphFile, SizeType vertices = 0, bool periodic = false)
	    : graphFile_(graphFile), vertices_(vertices)
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
		stripTrailingChars(data);
		loadGraphFromFile(data);
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

private:

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
	void loadGraphFromFile(PsimagLite::String data)
	{
		PsimagLite::String str;
		SizeType ind = readUntil(str, 0, data, '\n');

		VectorStringType tokens;
		PsimagLite::split(tokens, str, ":");
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
	VectorVectorSizeType allNeighbors_;
};
}
#endif // GRAPH_HH
