#ifndef _WEB_UI_HPP_
#define _WEB_UI_HPP_

#include "WebUIOutput.hpp"

// WUI outputs a to BOTH cout and the app...
#define WUI(KIND,MESSAGE) {                                           \
    WebUI::webUIOutput.message_flag = true;                           \
    if (mpi_rank==0) {                                                \
      std::cout << MESSAGE << std::endl;                              \
      stringstream _ss_;                                                \
      _ss_ << MESSAGE;                                                  \
      WebUI::webUIOutput.add_(UIMessage_(KIND,_ss_.str()));             \
    }                                                                 \
  }

namespace WebUI {

  extern WebUIOutput webUIOutput;
  
}

#endif
