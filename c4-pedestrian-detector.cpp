#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <string>
#include <cmath>
#include <vector>

#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#define USE_DOUBLE

#ifdef USE_DOUBLE
typedef double REAL;
#else
typedef float REAL;
#endif

template<class T> class Array2dC;

template<class T>
class Array2d
{
public:
    int nrow;
    int ncol;
    T** p;
public:
    Array2d():nrow(0),ncol(0),p(NULL) { }
    Array2d(const int nrow,const int ncol):nrow(0),ncol(0),p(NULL)
    {
        Create(nrow,ncol);
    }
    Array2d(const Array2d<T>& source);
    virtual ~Array2d()
    {
        Clear();
    }

    Array2d<T>& operator=(const Array2d<T>& source);
    void Create(const int _nrow,const int _ncol);
    void Swap(Array2d<T>& array2);
    void Clear();
    void Zero(const T t = 0);
};

template<class T>
class Array2dC
{
public:
    int nrow;
    int ncol;
    T** p;
    T* buf;
public:
    Array2dC():nrow(0),ncol(0),p(NULL),buf(NULL) {}
    Array2dC(const int nrow,const int ncol):nrow(0),ncol(0),p(NULL),buf(NULL)
    {
        Create(nrow,ncol);
    }
    Array2dC(const Array2dC<T>& source);
    virtual ~Array2dC()
    {
        Clear();
    }

    Array2dC<T>& operator=(const Array2dC<T>& source);
    void Create(const int _nrow,const int _ncol);
    void Swap(Array2dC<T>& array2);
    void Zero(const T t = 0);
    void Clear();
};

template<class T>
Array2d<T>::Array2d(const Array2d<T>& source):nrow(0),ncol(0),p(NULL)
{
    if(source.p!=NULL)
    {
        Create(source.nrow,source.ncol);
        for(int i=0; i<nrow; i++) std::copy(source.p[i],source.p[i]+ncol,p[i]);
    }
}

template<class T>
Array2d<T>& Array2d<T>::operator=(const Array2d<T>& source)
{
    if(source.p!=NULL)
    {
        Create(source.nrow,source.ncol);
        for(int i=0; i<nrow; i++) std::copy(source.p[i],source.p[i]+ncol,p[i]);
    }
    else
        Clear();
    return *this;
}

template<class T>
void Array2d<T>::Create(const int _nrow,const int _ncol)
{
    assert(_nrow>0 && _ncol>0);
    Clear();
    nrow = _nrow;
    ncol = _ncol;
    p = new T*[nrow];
    assert(p!=NULL);
    for(int i=0; i<nrow; i++)
    {
        p[i] = new T[ncol];
        assert(p[i]!=NULL);
    }
}

template<class T>
void Array2d<T>::Swap(Array2d<T>& array2)
{
    std::swap(nrow,array2.nrow);
    std::swap(ncol,array2.ncol);
    std::swap(p,array2.p);
}

template<class T>
void Array2d<T>::Zero(const T t)
{
    if(nrow>0)
    {
        for(int i=0; i<nrow; i++) std::fill(p[i],p[i]+ncol,t);
    }
}

template<class T>
void Array2d<T>::Clear()
{
    for(int i=0; i<nrow; i++)
    {
        delete[] p[i];
        p[i] = NULL;
    }
    delete[] p;
    p = NULL;
    nrow = ncol = 0;
}

template<class T>
Array2dC<T>::Array2dC(const Array2dC<T>& source):nrow(0),ncol(0),p(NULL),buf(NULL)
{
    if(source.buf!=NULL)
    {
        Create(source.nrow,source.ncol);
        std::copy(source.buf,source.buf+nrow*ncol,buf);
    }
}

template<class T>
Array2dC<T>& Array2dC<T>::operator=(const Array2dC<T>& source)
{
    if(source.buf!=NULL)
    {
        Create(source.nrow,source.ncol);
        std::copy(source.buf,source.buf+nrow*ncol,buf);
    }
    else
        Clear();
    return *this;
}

template<class T>
void Array2dC<T>::Create(const int _nrow,const int _ncol)
{
    assert(_nrow>0 && _ncol>0);
    if(nrow==_nrow && ncol==_ncol) return;
    Clear();
    nrow = _nrow;
    ncol = _ncol;
    buf = new T[nrow*ncol];
    assert(buf!=NULL);
    p = new T*[nrow];
    assert(p!=NULL);
    for(int i=0; i<nrow; i++) p[i] = buf + i * ncol;
}

template<class T>
void Array2dC<T>::Swap(Array2dC<T>& array2)
{
    std::swap(nrow,array2.nrow);
    std::swap(ncol,array2.ncol);
    std::swap(p,array2.p);
    std::swap(buf,array2.buf);
}

template<class T>
void Array2dC<T>::Zero(const T t)
{
    if(nrow>0) std::fill(buf,buf+nrow*ncol,t);
}

template<class T>
void Array2dC<T>::Clear()
{
    delete[] buf;
    buf = NULL;
    delete[] p;
    p = NULL;
    nrow = ncol = 0;
}


/*****************************************/
// IntImage.h
/*****************************************/

template<class T>
class IntImage:public Array2dC<T>
{
private:
    IntImage(const IntImage<T> &source) { } // prohibit copy constructor

public:
    IntImage():variance(0.0),label(-1) { }
    virtual ~IntImage()
    {
        Clear();
    }

    virtual void Clear(void);
    inline void SetSize(const int h, const int w);
    bool Load(cv::Mat img, const char channel = 'I');
    void Save(const std::string& filename) const;
    void Swap(IntImage<T>& image2);

    void CalcIntegralImageInPlace(void);
    void Resize(IntImage<T> &result,const REAL ratio) const;
    void Resize(IntImage<T>& result,const int height,const int width) const;

    IntImage<T>& operator=(const IntImage<T>& source);

    void Sobel(IntImage<REAL>& result);
public:
    using Array2dC<T>::nrow;
    using Array2dC<T>::ncol;
    using Array2dC<T>::buf;
    using Array2dC<T>::p;
    REAL variance;
    int label;
};

template<class T>
void IntImage<T>::Clear(void)
{
    Array2dC<T>::Clear();
    variance = 0.0;
    label = -1;
}

template<class T>
bool IntImage<T>::Load(cv::Mat img, const char channel)
{
    if (img.empty()) return false;

    if (channel == 'R' || channel == 'G' || channel == 'B')
    {
        int c;
        if (channel == 'B') c = 0;
        else if (channel == 'G') c = 1;
        else c = 2; // OpenCV is 'BGR' ordering
        cv::Mat planes[3];
        split(img, planes);
        img = planes[c];
    }
    else // use gray scale for all others
    {
        cv::cvtColor(img, img, cv::COLOR_BGR2GRAY);
    }

    SetSize(img.rows, img.cols);
    for(int i=0,ih=img.rows,iw=img.cols; i<ih; i++)
    {
        T* pdata = p[i];
        unsigned char* pimg = reinterpret_cast<unsigned char*>(img.data+img.step*i);
        for(int j=0; j<iw; j++) pdata[j] = pimg[j];
    }
    if(Show_Detection_Steps)
        imshow("1-Load", img);
    return true;
}

template<class T>
void IntImage<T>::Save(const std::string& filename) const
{
    cv::Mat img(cv::Size(ncol,nrow),CV_8U);
    for(int r=0; r<img.rows; r++)
    {
        T* pdata = p[r];
        unsigned char* pimg = img.ptr<unsigned char>(r);
        for(int c=0; c<img.cols; c++)
            pimg[c] = (unsigned char)pdata[c];
    }
    cv::imwrite(filename.c_str(),img);
}

template<class T>
void IntImage<T>::SetSize(const int h,const int w)
{
    if((h == nrow) && (w == ncol)) return;
    Clear();
    Array2dC<T>::Create(h,w);
}

template<class T>
IntImage<T>& IntImage<T>::operator=(const IntImage<T>& source)
{
    if(&source==this) return *this;
    SetSize(source.nrow,source.ncol);
    std::copy(source.buf,source.buf+nrow*ncol,buf);
    label = source.label;
    variance = source.variance;
    return *this;
}

template<class T>
void IntImage<T>::Resize(IntImage<T> &result,const REAL ratio) const
{
    Resize(result,int(nrow*ratio),int(ncol*ratio));
}

template<class T>
void IntImage<T>::Resize(IntImage<T>& result,const int height,const int width) const
{
    assert(height>0 && width>0);
    result.SetSize(height,width);
    REAL ixratio = nrow*1.0/height, iyratio = ncol*1.0/width;

    REAL* p_y = new REAL[result.ncol];
    assert(p_y!=NULL);
    int* p_y0 = new int[result.ncol];
    assert(p_y0!=NULL);
    for(int i=0; i<width; i++)
    {
        p_y[i] = i*iyratio;
        p_y0[i] = (int)p_y[i];
        if(p_y0[i]==ncol-1) p_y0[i]--;
        p_y[i] -= p_y0[i];
    }

    for(int i=0; i<height; i++)
    {
        int x0;
        REAL x;
        x = i*ixratio;
        x0 = (int)x;
        if(x0==nrow-1) x0--;
        x -= x0;
        T* rp = result.p[i];
        const T* px0 = p[x0];
        const T* px1 = p[x0+1];
        for(int j=0; j<width; j++)
        {
            int y0=p_y0[j];
            REAL y=p_y[j],fx0,fx1;

            fx0 = REAL(px0[y0] + y*(px0[y0+1]-px0[y0]));
            fx1 = REAL(px1[y0] + y*(px1[y0+1]-px1[y0]));

            rp[j] = T(fx0 + x*(fx1-fx0));
        }
    }

    delete[] p_y;
    p_y=NULL;
    delete[] p_y0;
    p_y0=NULL;
}

template<class T>
void IntImage<T>::CalcIntegralImageInPlace(void)
// We pad a zero column and a zero row, so 24*24 image will be 25*25 in size
// if the input image is not padded, the results on 1st row will be problematic
{
    for(int i=1; i<ncol; i++)   // process the first line
    {
        buf[i] += buf[i-1];
    }
    for(int i=1; i<nrow; i++)
    {
        REAL partialsum = 0;
        T* curp = p[i];
        T* prep = p[i-1];
        for(int j=0; j<ncol; j++)
        {
            partialsum += REAL(curp[j]);
            curp[j] = prep[j] + partialsum;
        }
    }
}

template<class T>
void IntImage<T>::Swap(IntImage<T>& image2)
{
    Array2dC<T>::Swap(image2);
    std::swap(variance,image2.variance);
    std::swap(label,image2.label);
}

template<class T>
void IntImage<T>::Sobel(IntImage<REAL>& result)
{
    // compute the Sobel gradient. For now, we just use the very inefficient way. Optimization can be done later
    for(int i=1; i<nrow-1; i++)
    {
        T* p1 = p[i-1];
        T* p2 = p[i];
        T* p3 = p[i+1];
        REAL* pr = result.p[i];
        for(int j=1; j<ncol-1; j++)
        {
            REAL gx =     p1[j-1] - p1[j+1]
                          + 2*(p2[j-1]   - p2[j+1])
                          +    p3[j-1] - p3[j+1];
            REAL gy =     p1[j-1] - p3[j-1]
                          + 2*(p1[j]   - p3[j])
                          +    p1[j+1] - p3[j+1];
            pr[j] = gx*gx+gy*gy;
        }
    }
}


/*****************************************/
// Pedestrian.h
/*****************************************/

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
        for(int i=0; i<length; i++) delete nodes[i];
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

/*****************************************/
// Pedestrian_ICRA.cpp
/*****************************************/

const int HUMAN_height = 108;
const int HUMAN_width = 36;
const int HUMAN_xdiv = 9;
const int HUMAN_ydiv = 4;
static const int EXT = 1;

// The detector
DetectionScanner scanner(HUMAN_height,HUMAN_width,HUMAN_xdiv,HUMAN_ydiv,256,0.7);
bool Show_Detection_Steps;
// ---------------------------------------------------------------------
// Helper functions


// compute the Sobel image "ct" from "original"
void ComputeCT(IntImage<double>& original,IntImage<int>& ct)
{
    ct.Create(original.nrow,original.ncol);
    for(int i=2; i<original.nrow-2; i++)
    {
        double* p1 = original.p[i-1];
        double* p2 = original.p[i];
        double* p3 = original.p[i+1];
        int* ctp = ct.p[i];
        for(int j=2; j<original.ncol-2; j++)
        {
            int index = 0;
            if(p2[j]<=p1[j-1]) index += 0x80;
            if(p2[j]<=p1[j]) index += 0x40;
            if(p2[j]<=p1[j+1]) index += 0x20;
            if(p2[j]<=p2[j-1]) index += 0x10;
            if(p2[j]<=p2[j+1]) index += 0x08;
            if(p2[j]<=p3[j-1]) index += 0x04;
            if(p2[j]<=p3[j]) index += 0x02;
            if(p2[j]<=p3[j+1]) index ++;
            ctp[j] = index;
        }
    }
}

// Load SVM models -- linear SVM trained using LIBLINEAR
double UseSVM_CD_FastEvaluationStructure(const char* modelfile,const int m,Array2dC<double>& result)
{
    std::ifstream in(modelfile);
    if(in.good()==false)
    {
        std::cout<<"SVM model "<<modelfile<<" can not be loaded."<<std::endl;
        exit(-1);
    }
    std::string buffer;
    std::getline(in,buffer); // first line
    std::getline(in,buffer); // second line
    std::getline(in,buffer); // third line
    in>>buffer;
    assert(buffer=="nr_feature");
    int num_dim;
    in>>num_dim;
    assert(num_dim>0 && num_dim==m);
    std::getline(in,buffer); // end of line 4
    in>>buffer;
    assert(buffer=="bias");
    int bias;
    in>>bias;
    std::getline(in,buffer); //end of line 5;
    in>>buffer;
    assert(buffer=="w");
    std::getline(in,buffer); //end of line 6
    result.Create(1,num_dim);
    for(int i=0; i<num_dim; i++) in>>result.buf[i];
    double rho = 0;
    if(bias>=0) in>>rho;
    in.close();
    return rho;
}

// Load SVM models -- Histogram Intersectin Kernel SVM trained by libHIK
double UseSVM_CD_FastEvaluationStructure(const char* modelfile, const int m, const int upper_bound, Array2dC<double>& result)
{

    std::ifstream fs(modelfile, std::fstream::binary);
	if( !fs.is_open() )
	{
		std::cout << "SVM model " << modelfile << " can not be loaded." << std::endl;
		exit(-1);
	}
    // Header
    int rows, cols, type, channels;
    fs.read((char*)&rows, sizeof(int));         // rows
    fs.read((char*)&cols, sizeof(int));         // cols
    fs.read((char*)&type, sizeof(int));         // type
    fs.read((char*)&channels, sizeof(int));     // channels

    // Data
    cv::Mat mat(rows, cols, type);
    fs.read((char*)mat.data, CV_ELEM_SIZE(type) * rows * cols);

    int num_dim = m;

    result.Create(num_dim, upper_bound);
    for(int i=0; i<num_dim; i++)
        for (int j = 0; j < upper_bound; j++)
        {
            result.p[i][j]= mat.at<double>(i, j);
        }

    return -0.00455891;
}

// find the intersection of "this" and "rect2", and put into "result"
bool CRect::Intersect(CRect& result,const CRect& rect2) const
{
    if( Empty() || rect2.Empty() ||
            left >= rect2.right || rect2.left >= right ||
            top >= rect2.bottom || rect2.top >= bottom )
    {
        result.Clear();
        return false;
    }
    result.left   = std::max( left, rect2.left );
    result.right  = std::min( right, rect2.right );
    result.top    = std::max( top, rect2.top );
    result.bottom = std::min( bottom, rect2.bottom );
    return true;
}

// find the union of "this" and "rect2", and put into "result"
bool CRect::Union(CRect& result,const CRect& rect2) const
{
    if(Empty())
    {
        if(rect2.Empty())
        {
            result.Clear();
            return false;
        }
        else
            result = rect2;
    }
    else
    {
        if(rect2.Empty())
            result = *this;
        else
        {
            result.left   = std::min( left, rect2.left );
            result.right  = std::max( right, rect2.right );
            result.top    = std::min( top, rect2.top );
            result.bottom = std::max( bottom, rect2.bottom );
        }
    }
    return true;
}

// A simple post-process (NMS, non-maximal suppression)
// "result" -- rectangles before merging
//          -- after this function it contains rectangles after NMS
// "combine_min" -- threshold of how many detection are needed to survive
void PostProcess(std::vector<CRect>& result,const int combine_min)
{
    std::vector<CRect> res1;
    std::vector<CRect> resmax;
    std::vector<int> res2;
    bool yet;
    CRect rectInter;

    for(size_t i=0,size_i=result.size(); i<size_i; i++)
    {
        yet = false;
        CRect& result_i = result[i];
        for(size_t j=0,size_r=res1.size(); j<size_r; j++)
        {
            CRect& resmax_j = resmax[j];
            if(result_i.Intersect(rectInter,resmax_j))
            {
                if(  rectInter.Size()>0.6*result_i.Size()
                        && rectInter.Size()>0.6*resmax_j.Size()
                  )
                {
                    CRect& res1_j = res1[j];
                    resmax_j.Union(resmax_j,result_i);
                    res1_j.bottom += result_i.bottom;
                    res1_j.top += result_i.top;
                    res1_j.left += result_i.left;
                    res1_j.right += result_i.right;
                    res2[j]++;
                    yet = true;
                    break;
                }
            }
        }
        if(yet==false)
        {
            res1.push_back(result_i);
            resmax.push_back(result_i);
            res2.push_back(1);
        }
    }

    for(size_t i=0,size=res1.size(); i<size; i++)
    {
        const int count = res2[i];
        CRect& res1_i = res1[i];
        res1_i.top /= count;
        res1_i.bottom /= count;
        res1_i.left /= count;
        res1_i.right /= count;
    }

    result.clear();
    for(size_t i=0,size=res1.size(); i<size; i++)
        if(res2[i]>combine_min)
            result.push_back(res1[i]);
}

// If one detection (after NMS) is inside another, remove the inside one
void RemoveCoveredRectangles(std::vector<CRect>& result)
{
    std::vector<bool> covered;
    covered.resize(result.size());
    std::fill(covered.begin(),covered.end(),false);
    CRect inter;
    for(unsigned int i=0; i<result.size(); i++)
    {
        for(unsigned int j=i+1; j<result.size(); j++)
        {
            result[i].Intersect(inter,result[j]);
            double isize = inter.Size();
            if(isize>result[i].Size()*0.65)
                covered[i] = true;
            if(isize>result[j].Size()*0.65)
                covered[j] = true;
        }
    }
    std::vector<CRect> newresult;
    for(unsigned int i=0; i<result.size(); i++)
        if(covered[i]==false)
            newresult.push_back(result[i]);
    result.clear();
    result.insert(result.begin(),newresult.begin(),newresult.end());
    newresult.clear();
}

// End of Helper functions
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Functions that load the two classifiers
void LoadCascade(DetectionScanner& ds)
{
    std::vector<NodeDetector::NodeType> types;
    std::vector<int> upper_bounds;
    std::vector<std::string> filenames;

    types.push_back(NodeDetector::CD_LIN); // first node
    upper_bounds.push_back(100);
    filenames.push_back("combined.txt.model");
    types.push_back(NodeDetector::CD_HIK); // second node
    upper_bounds.push_back(353);
    filenames.push_back("combined.txt.model_");

    ds.LoadDetector(types,upper_bounds,filenames);
    // You can adjust these parameters for different speed, accuracy etc
    ds.cascade->nodes[0]->thresh += 0.8;
    ds.cascade->nodes[1]->thresh -= 0.095;
}

void DetectionScanner::LoadDetector(std::vector<NodeDetector::NodeType>& types,std::vector<int>& upper_bounds,std::vector<std::string>& filenames)
{
    size_t depth = types.size();
    assert(depth>0 && depth==upper_bounds.size() && depth==filenames.size());
    if(cascade)
        delete cascade;
    cascade = new CascadeDetector;
    assert(xdiv>0 && ydiv>0);
    for(unsigned int i=0; i<depth; i++)
        cascade->AddNode(types[i],(xdiv-EXT)*(ydiv-EXT)*baseflength,upper_bounds[i],filenames[i].c_str());

    hist.Create(1,baseflength*(xdiv-EXT)*(ydiv-EXT));
}

void NodeDetector::Load(const NodeType _type,const int _featurelength,const int _upper_bound,const int _index,const char* _filename)
{
    type = _type;
    index = _index;
    filename = _filename;
    featurelength = _featurelength;
    upper_bound = _upper_bound;
    if(type==CD_LIN)
        thresh = UseSVM_CD_FastEvaluationStructure(_filename,_featurelength,classifier);
    else if(type==CD_HIK)
        thresh = UseSVM_CD_FastEvaluationStructure(_filename,_featurelength,upper_bound,classifier);

    if(type==CD_LIN) type = LINEAR;
    if(type==CD_HIK) type = HISTOGRAM;
}

void CascadeDetector::AddNode(const NodeDetector::NodeType _type,const int _featurelength,const int _upper_bound,const char* _filename)
{
    if(length==size)
    {
        int newsize = size * 2;
        NodeDetector** p = new NodeDetector*[newsize];
        assert(p!=NULL);
        std::copy(nodes,nodes+size,p);
        size = newsize;
        delete[] nodes;
        nodes = p;
    }
    nodes[length] = new NodeDetector(_type,_featurelength,_upper_bound,length,_filename);
    length++;
}

// End of functions that load the two classifiers
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Detection functions

// initialization -- compute the Census Tranform image for CENTRIST
void DetectionScanner::InitImage(IntImage<double>& original)
{
    image = original;
    sobel.Create(image.nrow, image.ncol);
    image.Sobel(sobel);
    ComputeCT(sobel,ct);

    if (Show_Detection_Steps)
    {
        cv::Mat img_sobel(sobel.nrow, sobel.ncol, CV_64F, sobel.buf),img;
        img_sobel.convertTo(img, CV_8U, 1. / 256);
        imshow("2-Sobel", img);
        cv::Mat img_ct(ct.nrow, ct.ncol, CV_32S, ct.buf);
        img_ct.convertTo(img, CV_8U);
        imshow("3-ct", img);
    }

}

// combine the (xdiv-1)*(ydiv-1) integral images into a single one
void DetectionScanner::InitIntegralImages(const int stepsize)
{
    if(cascade->nodes[0]->type!=NodeDetector::LINEAR)
        return; // No need to prepare integral images

    const int hd = height/xdiv*2-2;
    const int wd = width/ydiv*2-2;
    scores.Create(ct.nrow,ct.ncol);
    scores.Zero(cascade->nodes[0]->thresh/hd/wd);
    double* linearweights = cascade->nodes[0]->classifier.buf;
    for(int i=0; i<xdiv-EXT; i++)
    {
        const int xoffset = height/xdiv*i;
        for(int j=0; j<ydiv-EXT; j++)
        {
            const int yoffset = width/ydiv*j;
            for(int x=2; x<ct.nrow-2-xoffset; x++)
            {
                int* ctp = ct.p[x+xoffset]+yoffset;
                double* tempp = scores.p[x];
                for(int y=2; y<ct.ncol-2-yoffset; y++)
                    tempp[y] += linearweights[ctp[y]];
            }
            linearweights += baseflength;
        }
    }
    scores.CalcIntegralImageInPlace();
    for(int i=2; i<ct.nrow-2-height; i+=stepsize)
    {
        double* p1 = scores.p[i];
        double* p2 = scores.p[i+hd];
        for(int j=2; j<ct.ncol-2-width; j+=stepsize)
            p1[j] += (p2[j+wd] - p2[j] - p1[j+wd]);
    }

    if (Show_Detection_Steps)
    {
        cv::Mat img_scores(scores.nrow, scores.ncol, CV_64F, scores.buf), img,cimg;
        cv::Mat img_sobel(sobel.nrow, sobel.ncol, CV_64F, sobel.buf);
        img_sobel.convertTo(img, CV_8U, 1. / 256);
        cv::cvtColor(img, cimg, cv::COLOR_GRAY2BGR);
        img_scores.convertTo(img, CV_8U, 187);
        cv::insertChannel(img, cimg, 2);
        imshow("2-Sobel", cimg);
        imshow("4-scores", img);

        int key = cv::waitKey();

        if (key == 's' || key == 'S')
        {
            Show_Detection_Steps = !Show_Detection_Steps;
            cv::destroyAllWindows();
        }
    }
}

// Resize the input image and then re-compute Sobel image etc
void DetectionScanner::ResizeImage()
{
    image.Resize(sobel,ratio);
    image.Swap(sobel);
    image.Sobel(sobel);
    ComputeCT(sobel,ct);
}

// The function that does the real detection
int DetectionScanner::FastScan(IntImage<double>& original,std::vector<CRect>& results,const int stepsize)
{
    if(original.nrow<height+5 || original.ncol<width+5) return 0;
    const int hd = height/xdiv;
    const int wd = width/ydiv;
    InitImage(original);
    results.clear();

    hist.Create(1,baseflength*(xdiv-EXT)*(ydiv-EXT));

    NodeDetector* node = cascade->nodes[1];
    double** pc = node->classifier.p;
    int oheight = original.nrow, owidth = original.ncol;
    CRect rect;
    while(image.nrow>=height && image.ncol>=width)
    {
        InitIntegralImages(stepsize);
        for(int i=2; i+height<image.nrow-2; i+=stepsize)
        {
            const double* sp = scores.p[i];
            for(int j=2; j+width<image.ncol-2; j+=stepsize)
            {
                if(sp[j]<=0) continue;
                int* p = hist.buf;
                hist.Zero();
                for(int k=0; k<xdiv-EXT; k++)
                {
                    for(int t=0; t<ydiv-EXT; t++)
                    {
                        for(int x=i+k*hd+1; x<i+(k+1+EXT)*hd-1; x++)
                        {
                            int* ctp = ct.p[x];
                            for(int y=j+t*wd+1; y<j+(t+1+EXT)*wd-1; y++)
                                p[ctp[y]]++;
                        }
                        p += baseflength;
                    }
                }
                double score = node->thresh;
                for(int k=0; k<node->classifier.nrow; k++) score += pc[k][hist.buf[k]];
                if(score>0)
                {
                    rect.top = i*oheight/image.nrow;
                    rect.bottom = (i+height)*oheight/image.nrow;
                    rect.left = j*owidth/image.ncol;
                    rect.right = (j+width)*owidth/image.ncol;
                    results.push_back(rect);
                }
            }
        }
        ResizeImage();
    }
    return 0;
}

int main(int argc, char* argv[])
{
    std::cout << "usage:" << std::endl;
    std::cout << argv[0] << " <video_file>" << std::endl << std::endl;
    std::cout << "keys:" << std::endl;
    std::cout << "space : toggle using simple post-process (NMS, non-maximal suppression)" << std::endl;
    std::cout << "0     : waits to process next frame until a key pressed" << std::endl;
    std::cout << "1     : doesn't wait to process next frame" << std::endl;
    std::cout << "2     : resize frames 1/2" << std::endl;
    std::cout << "3     : don't resize frames" << std::endl;
    std::cout << "4     : resize frames 1/4" << std::endl;
    std::cout << "s     : toggle Show_Detection_Steps" << std::endl;
    std::cout << "r     : increase resize ratio" << std::endl;
    std::cout << "e     : decrease resize ratio" << std::endl;

    if (argc < 2)
        return 0;

    cv::TickMeter tm;
    cv::Animation animation;
    cv::Mat src, img;
    bool is_video = true;

    cv::VideoCapture capture(argv[1]);
    if (capture.get(cv::CAP_PROP_FRAME_COUNT) == 1)
    {
        capture >> img;
        is_video = false;
    }


    LoadCascade(scanner);
    std::cout << "Detectors loaded." << std::endl;
    int key = 0;
    int wait_time = 1;
    float fx = 1;
    bool rect_organization = true;
    IntImage<double> original;

    while (key != 27)
    {
        if (is_video)
            capture >> src;
        else
            src = img.clone();

        if (src.empty()) break;

        if (fx < 1)
        {
            cv::resize(src, src, cv::Size(), fx, fx);
        }

        original.Load(src);
        std::vector<CRect> results;
        tm.reset();
        tm.start();
        scanner.FastScan(original, results, 2);

        if (rect_organization)
        {
            PostProcess(results, 2);
            PostProcess(results, 0);
            RemoveCoveredRectangles(results);
        }
        tm.stop();
        std::cout << "\ndetection time :" << tm.getTimeMilli() << " ms.";

        for (size_t i = 0; i < results.size(); i++)
        {
            cv::rectangle(src, cv::Point2d(results[i].left, results[i].top), cv::Point2d(results[i].right, results[i].bottom), cv::Scalar(0, 255, 0), 2);
        }

        cv::imwrite("result.jpg", src);
        cv::imshow("result", src);
        cv::Mat resized;
        cv::resize(src, resized, cv::Size(), 0.5, 0.5);
        animation.frames.push_back(resized);
        animation.durations.push_back(cvRound(tm.getTimeMilli()/2));
        key = cv::waitKey(wait_time);

        if (key == 32)
            rect_organization = !rect_organization;

        if (key == 48)
            wait_time = 0;

        if (key == 49)
            wait_time = 1;

        if (key == 50)
            fx = 0.5;

        if (key == 51)
            fx = 1;

        if (key == 52)
            fx = 0.25;

        if (key == 's' || key == 'S')
        {
            Show_Detection_Steps = !Show_Detection_Steps;
            cv::destroyAllWindows();
        }

        if (key == 'r' || key == 'R')
        {
            if (scanner.ratio < 0.9)
                scanner.ratio *= 1.1;
            else scanner.ratio = 0.95;
            std::cout << "scanner.ratio :" << scanner.ratio << std::endl;
        }

        if (key == 'e' || key == 'E')
        {
            if (scanner.ratio > 0.2)
                scanner.ratio *= 0.9;
            std::cout << "scanner.ratio :" << scanner.ratio << std::endl;
        }

    }
    tm.reset();
    tm.start();
    cv::imwriteanimation("result.png", animation);
    tm.stop();
    std::cout << "\nsave time :" << tm.getTimeMilli() << " ms.";
    tm.reset();
    tm.start();
    cv::imwriteanimation("result_IMWRITE_PNG_COMPRESSION_6.png", animation,{cv::IMWRITE_PNG_COMPRESSION,6});
    tm.stop();
    std::cout << "\nsave time :" << tm.getTimeMilli() << " ms.";
    tm.reset();
    tm.start();
    cv::imwriteanimation("result.webp", animation);
    tm.stop();
    std::cout << "\nsave time :" << tm.getTimeMilli() << " ms.";
    tm.reset();
    tm.start();
    cv::imwriteanimation("result.avif", animation);
    tm.stop();
    std::cout << "\nsave time :" << tm.getTimeMilli() << " ms.";
    tm.reset();
    tm.start();
    cv::imwriteanimation("result.gif", animation);
    tm.stop();
    std::cout << "\nsave time :" << tm.getTimeMilli() << " ms.";
    return 0;
}
