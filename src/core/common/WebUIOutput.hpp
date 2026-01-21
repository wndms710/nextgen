#ifndef _WEB_UI_OUTPUT_HPP_
#define _WEB_UI_OUTPUT_HPP_

#include "CTI.hpp"

using namespace CTI;

enum ui_message_t {
  MESSAGE,
  INFO,
  WARN,
  ERR,
  STITCH_IN,
  CHARLES_IN
};

class UIMessage_ {
private:
  ui_message_t type;
  bool closable;  // user closes notification
  int timeout;  // notification dissappears (specified in ms)
  vector<string> lines;

public:
  UIMessage_() {
    type = MESSAGE;
    closable = false;
    timeout = 5000;
  }

  UIMessage_(const ui_message_t _type) {
    setType(_type);
    closable = false;
    timeout = 5000;
  }

  /*
  UIMessage(const ui_message_t _type,const string line) {
    setType(_type);
    closable = false;
    timeout = 5000;
    lines.push_back(line);
  }
  */

  UIMessage_(const ui_message_t _type,const string message) {
    setType(_type);
    closable = false;
    timeout = 5000;
    std::size_t pos = message.find_first_of('\n',0);
    if (pos == std::string::npos) {
      // no returns, so just add amessage as a line...
      lines.push_back(message);
    }
    else {
      // add the first part up to the return...
      lines.push_back(message.substr(0,pos));
      std::size_t pos2 = message.find_first_of('\n',pos+1);
      while (pos2 != std::string::npos) {
        // keep adding up to the next returns...
        lines.push_back(message.substr(pos+1,pos2-pos-1));
        pos = pos2;
        pos2 = message.find_first_of('\n',pos+1);
      }
      // if the line doesn't end in return, there is one more to add...
      if (pos+1 < message.length())
        lines.push_back(message.substr(pos+1));
    }
  }

  UIMessage_(const ui_message_t _type,const vector<string> _lines) {
    setType(_type);
    closable = false;
    timeout = 5000;
    lines = _lines;
  }

  UIMessage_(const UIMessage_& other) {
    type = other.type;
    closable = other.closable;
    timeout = other.timeout;
    lines = other.lines;
  }

  ~UIMessage_() {}

  inline int size() const {
    return lines.size();
  }
  inline string& operator[] (const int index) {
    return lines[index];
  }

  void clear() {
    lines.clear();
  }

  void addLine(const string line) {
    lines.push_back(line);
  }
  void addLine(ostream& ss) {
    lines.push_back(dynamic_cast<stringstream&>(ss).str());
  }

  void dumpLines() {
    for (vector<string>::iterator it=lines.begin(); it!=lines.end(); ++it) {
      COUT1(*it);
    }
  }

  void setType(const ui_message_t _type) {
    if (_type >= 0 && _type <= 5) type = _type;
  }
  inline ui_message_t& getType() {
    return type;
  }

  void setClosable(const bool _closable) {
    closable = _closable;
  }
  inline bool& isClosable() {
    return closable;
  }

  void setTimeout(const int ms) {
    timeout = ms;
  }
  inline int& getTimeout() {
    return timeout;
  }
};

class WebUIOutput {
private:
  vector<UIMessage_> messages;
  //bool writable;
  bool rebuildMenu;

public:
  bool newImage; // TODO : get rid of this, and determine if an image is required from the bool of all the messages
  bool message_flag;
  WebUIOutput() {
    newImage = false; // default behavior
    message_flag = false;
    rebuildMenu = false;
    //writable = true;
  }

  /*
  void writeOK(const bool _writable) {
    writable = _writable;
  }
  */
  void add_(const UIMessage_ message) {
    //if (writable) messages.push_back(message);
    message_flag = true;
    messages.push_back(message);
  }

  inline UIMessage_& operator[] (const int index) {
    if (index < 0 || index >= int(messages.size())) CWARN("out-of-bounds access into webUIOutput["<<index<<"]");
    return messages[index];
  }
  inline int size() const {
    return messages.size();
  }
  void ensureImage() {
    newImage = true;
  }
  bool skipImage() {
    return !newImage;
  }
  void ensureMenu() {
    rebuildMenu = true;
  }
  bool buildMenu() {
    return rebuildMenu;
  }
  void clear() {
    messages.clear();
    // AND reset default behavior for next JSON write...
    newImage = false;
    rebuildMenu = false;
    message_flag = false;
  }
  bool empty() const {
    return messages.empty();
  }
  void writeJson(FILE * fp) {
    if (!empty()) {
      fprintf(fp,",\n\"webUIOutput\": [");
      bool first_time = true;
      for (int i=0, limit=messages.size(); i<limit; ++i) {
        if (!first_time) fprintf(fp,",");
        fprintf(fp,"\n  {");  // open message
        // specify type
        switch (messages[i].getType()) {
        case 0:
          fprintf(fp,"\n    \"type\": \"message\",");
          break;
        case 1:
          fprintf(fp,"\n    \"type\": \"info\",");
          break;
        case 2:
          fprintf(fp,"\n    \"type\": \"warning\",");
          break;
        case 3:
          fprintf(fp,"\n    \"type\": \"error\",");
          break;
        case 4:
          fprintf(fp,"\n    \"type\": \"s_input_file\",");
          break;
        case 5:
          fprintf(fp,"\n    \"type\": \"c_input_file\",");
          break;
        }
        fprintf(fp,"\n    \"closable\": %s,",((messages[i].isClosable()) ? "true":"false"));
        fprintf(fp,"\n    \"timeout\": %d,",messages[i].getTimeout());
        // message as a list of strings
        fprintf(fp,"\n    \"message\": [");
        bool first_line = true;
        for (int line=0, nlines=messages[i].size(); line<nlines; ++line) {
          if (!first_line) fprintf(fp,",");
          fprintf(fp,"\n      \"%s\"",messages[i][line].c_str());
          if (first_line) first_line = false;
        }
        fprintf(fp,"\n    ]");  // end lines
        fprintf(fp,"\n  }");  // close message
        if (first_time) first_time = false;
      }
      fprintf(fp,"\n]");
    }
    if (!newImage) {
      COUT1(" > alerting App to skip image");
      fprintf(fp,",\n\"newImage\":false");
    }
    if (rebuildMenu) {
      COUT1(" > alerting App to rebuild menus");
      fprintf(fp,",\n\"rebuildMenu\":true");
    }
  }

};

#endif
