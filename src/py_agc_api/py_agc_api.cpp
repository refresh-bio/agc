#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "../core/agc_decompressor.h"
#include "../lib-cxx/agc-api.h"

//binding STL container std::vector<std::string>
PYBIND11_MAKE_OPAQUE(std::vector<std::string>);

namespace py = pybind11;


PYBIND11_MODULE(py_agc_api, m) {
	m.doc() = "Python wrapper for AGC_API."; // optional module docstring
	
    // StringVector can be used in Python code (binding std::vector<std::string>)
    py::bind_vector<std::vector<std::string>>(m, "StringVector")
        .def(py::init<>())
        .def("clear", &std::vector<std::string>::clear)
        .def("pop_back", &std::vector<std::string>::pop_back)
        .def("__len__", [](const std::vector<std::string> &v) { return v.size(); })
        .def("__iter__", [](std::vector<std::string> &v) {
               return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>());
    
    // Class that represents agc archive
    py::class_<CAGCFile>(m, "CAGCFile")
        .def(py::init<>()) //parameterless constructor
        
        //Open(file_name, prefetching = true) opens agc archive
        //
        //@return true for success and false for error
        .def("Open", &CAGCFile::Open)
        
        //Close() closes opened archive
        //@return true for success and false for error
        .def("Close", &CAGCFile::Close) //Close() closes opened archive
        
        //NSample()
        //@returns  number of samples in the archive
        .def("NSample", &CAGCFile::NSample)
    
        //GetReferenceSample()
        //@returns reference sample
        .def("GetReferenceSample", [](CAGCFile& ptr){ std::string s;  ptr.GetReferenceSample(s); return s;})
             
    
        //NCtg(sample)
        //@returns  number of contig in sample
        .def("NCtg", &CAGCFile::NCtg)
    
        //ListSample(samples: StringVector)
        //@param samples  vector of strings (StringVector) with sample names (returned value)
        //@return number of samples to be written to
        .def("ListSample", &CAGCFile::ListSample)
    
        //ListCtg(sample, names: StringVector)
        //@param sample sample name
        //@param names  vector of strings (StringVector) with contig names (returned value)
        //@return number of contigs in the sample
        .def("ListCtg", &CAGCFile::ListCtg)
    
        //GetCtgLen(sample, name)
        //Get the length of a contig.
        //@param sample   sample name;
        //@param name     contig name
        //@return contig length, or <0 for errors
        .def("GetCtgLen", &CAGCFile::GetCtgLen)
        
        //GetCtgSeq(sample, name, start, end)
        //@param sample   sample name
        //@param name     contig name
        //@param start    start offset
        //@param end      end offset
        //@return contig sequence
        .def("GetCtgSeq", [](CAGCFile& ptr, const std::string& sample, const std::string& name, int start, int end) { std::string s;  ptr.GetCtgSeq(sample, name, start, end, s); return s;})
    
        //GetCtgSeq(name, start, end)
        //@param name     contig name
        //@param start    start offset
        //@param end      end offset
        //@return contig sequence (if unique name across all contigs in all samples)
        .def("GetCtgSeq", [](CAGCFile& ptr, const std::string& name, int start, int end) { std::string s; std::string empty; ptr.GetCtgSeq(empty, name, start, end, s); return s;})
    ;
		
}

