#ifndef __PESTRIAN_H_
#define __PESTRIAN_H_

#include <vector>
#include "IntImage.h"

// my replacement for the CRect class in MS MFC -- only provides a limited number of functions
class CRect
{
    public:
        double left;
        double top;
        double right;
        double bottom;
    public:
        CRect()
        {
            Clear();
        }
        ~CRect()
        {
            Clear();
        }
    public:
        bool Empty() const
        {
            return (left >= right) || (top >= bottom);
        }
        void Clear()
        {
            left = right = top = bottom = 0;
        }
        double Size() const
        {
            if(Empty())
                return 0;
            else
                return (bottom-top)*(right-left);
        }
    // Intersect and Union of two rectangles, both function should be able to run when &result==this
        bool Intersect(CRect& result,const CRect& rect2) const;
        bool Union(CRect& result,const CRect& rect2) const;
};

class NodeDetector
{
public:
	enum NodeType { CD_LIN, CD_HIK, LINEAR, HISTOGRAM };
public:
	int type; // linear or histogram?
	Array2dC<double> classifier;
	double thresh;
	int featurelength;
	int upper_bound;
	int index;
	std::string filename;
public:
	NodeDetector(const NodeType _type,const int _featurelength,const int _upper_bound,const int _index,const char* _filename)
	{
		Load(_type,_featurelength,_upper_bound,_index,_filename);
		minvalue = DBL_MAX;
		maxvalue = -minvalue;
	}
	~NodeDetector()
	{
	}
	
	void Load(const NodeType _type,const int _featurelength,const int _upper_bound,const int _index,const char* _filename);
	bool Classify(int* f);
private:
	double minvalue;
	double maxvalue;
public:
	void SetValues(const double v)
	{
		if(v>maxvalue) maxvalue = v;
		if(v<minvalue) minvalue = v;
	}
};

class CascadeDetector
{
public:
	int size;
	int length;
	NodeDetector** nodes;
public:
	
public:
	CascadeDetector()
		: size(20), length(0)
	{
		nodes = new NodeDetector*[size];
	}
	~CascadeDetector()
	{
		for(int i=0;i<length;i++) delete nodes[i];
		delete[] nodes;
	}
	
	void AddNode(const NodeDetector::NodeType _type,const int _featurelength,const int _upper_bound,const char* _filename);
};

class DetectionScanner // who does the dirty jobs
{
public:
	int height,width;
	int xdiv,ydiv;
	int baseflength;
	double ratio;
	
	CascadeDetector* cascade;
public:
	DetectionScanner()
		: height(0), width(0), xdiv(0), ydiv(0), baseflength(0), ratio(0.0), cascade(NULL), integrals(NULL)
	{
	}
	DetectionScanner(const int _height,const int _width,const int _xdiv,const int _ydiv,
					 const int _baseflength,const double _ratio)
		:height(_height),width(_width),xdiv(_xdiv),ydiv(_ydiv),
		 baseflength(_baseflength),ratio(_ratio),cascade(NULL),integrals(NULL)
	{
	}
	~DetectionScanner()
	{
		delete cascade;
		delete[] integrals;
	}
public:
	void LoadDetector(std::vector<NodeDetector::NodeType>& types,std::vector<int>& upper_bounds,std::vector<std::string>& filenames);
private:
	IntImage<double>* integrals;
	IntImage<double> image,sobel;
	IntImage<int> ct;
	Array2dC<int> hist;
	IntImage<double> scores;
	
	void InitImage(IntImage<double>& original);
	void InitIntegralImages(const int stepsize);
	void ResizeImage();
public:
	int Scan(IntImage<double>& original,std::vector<CRect>& results,const int stepsize,const int round,std::ofstream* out,const int upper_bound);
	int FastScan(IntImage<double>& original,std::vector<CRect>& results,const int stepsize);
	int FeatureLength() const
	{
		return (xdiv-1)*(ydiv-1)*baseflength;
	}
};

void RunFiles();

#endif // __PESTRIAN_H_
