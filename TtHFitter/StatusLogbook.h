#ifndef StatusLogbook_H_
#define StatusLogbook_H_

#include <string>

void WriteErrorStatus(const std::string&, const std::string&);
void WriteWarningStatus(const std::string&, const std::string&);
void WriteDebugStatus(const std::string&, const std::string&);
void WriteVerboseStatus(const std::string&, const std::string&);
void WriteInfoStatus(const std::string&, const std::string&);

#endif

