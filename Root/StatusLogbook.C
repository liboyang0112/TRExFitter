#include "TtHFitter/StatusLogbook.h"
#include "TtHFitter/Common.h"

#include <iostream>

int fMaxStrLen  = 54000;

void WriteErrorStatus(const std::string& classname, const std::string& info)
{

    int fNWhitespace = fMaxStrLen - 42 - info.size() - classname.size();
    
    std::string Whitespace = "";
    
    for(int i = 0; i < fNWhitespace; ++i)
        Whitespace += " ";
    
    std::string outputstring = "=== ERROR::"+classname+": "+info;
    
    // always print error
    std::cerr << "\033[1;31m" << outputstring.c_str() << "\33[0m" << std::endl;

}

void WriteWarningStatus(const std::string& classname, const std::string& info)
{

    int fNWhitespace = fMaxStrLen - 42 - info.size() - classname.size();
    
    std::string Whitespace = "";
    
    for(int i = 0; i < fNWhitespace; ++i)
        Whitespace += " ";
    
    std::string outputstring = "=== WARNING::"+classname+": "+info;
    
    // always print warnings
    std::cout << "\033[1;33m" << outputstring.c_str() << "\33[0m" << std::endl;

}

void WriteInfoStatus(const std::string& classname, const std::string& info)
{

    int fNWhitespace = fMaxStrLen - 42 - info.size() - classname.size();
    
    std::string Whitespace = "";
    
    for(int i = 0; i < fNWhitespace; ++i)
        Whitespace += " ";
    
    std::string outputstring = "=== INFO::"+classname+": "+info;
    
    if (TtHFitter::DEBUGLEVEL > 0) std::cout << "\033[1;32m" << outputstring.c_str() << "\33[0m" << std::endl;

}

void WriteDebugStatus(const std::string& classname, const std::string& info)
{

    int fNWhitespace = fMaxStrLen - 42 - info.size() - classname.size();
    
    std::string Whitespace = "";
    
    for(int i = 0; i < fNWhitespace; ++i)
        Whitespace += " ";
    
    std::string outputstring = "=== DEBUG::"+classname+": "+info;
    
    if (TtHFitter::DEBUGLEVEL > 1) std::cout << outputstring.c_str() << std::endl;

}
