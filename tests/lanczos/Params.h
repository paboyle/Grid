/*
  Params IO

  Author: Christoph Lehner
  Date:   2017
*/

#define PADD(p,X) p.get(#X,X);

class Params {
 protected:

  std::string trim(const std::string& sc) {
    std::string s = sc;
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
				    std::not1(std::ptr_fun<int, int>(std::isspace))));
    s.erase(std::find_if(s.rbegin(), s.rend(),
			 std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
  }

 public:

  std::map< std::string, std::string > lines;
  std::string _fn;

 Params(const char* fn) : _fn(fn) {
    FILE* f = fopen(fn,"rt");
    assert(f);
    while (!feof(f)) {
      char buf[4096];
      if (fgets(buf,sizeof(buf),f)) {
	if (buf[0] != '#' && buf[0] != '\r' && buf[0] != '\n') {
	  char* sep = strchr(buf,'=');
	  assert(sep);
	  *sep = '\0';
	  lines[trim(buf)] = trim(sep+1);
	}
      }
    }      
    fclose(f);
  }

  ~Params() {
  }

  std::string loghead() {
    return _fn + ": ";
  }

  bool has(const char* name) {
    auto f = lines.find(name);
    return (f != lines.end());
  }

  const std::string& get(const char* name) {
    auto f = lines.find(name);
    if (f == lines.end()) {
      std::cout << Grid::GridLogMessage << loghead() << "Could not find value for " << name << std::endl;
      abort();
    }
    return f->second;
  }

  void parse(std::string& s, const std::string& cval) {
    std::stringstream trimmer;
    trimmer << cval;
    s.clear();
    trimmer >> s;
  }

  void parse(int& i, const std::string& cval) {
    assert(sscanf(cval.c_str(),"%d",&i)==1);
  }

  void parse(long long& i, const std::string& cval) {
    assert(sscanf(cval.c_str(),"%lld",&i)==1);
  }

  void parse(double& f, const std::string& cval) {
    assert(sscanf(cval.c_str(),"%lf",&f)==1);
  }

  void parse(float& f, const std::string& cval) {
    assert(sscanf(cval.c_str(),"%f",&f)==1);
  }

  void parse(bool& b, const std::string& cval) {
    std::string lcval = cval;
    std::transform(lcval.begin(), lcval.end(), lcval.begin(), ::tolower);
    if (lcval == "true" || lcval == "yes") {
      b = true;
    } else if (lcval == "false" || lcval == "no") {
      b = false;
    } else {
      std::cout << "Invalid value for boolean: " << b << std::endl;
      assert(0);
    }
  }

  void parse(std::complex<double>& f, const std::string& cval) {
    double r,i;
    assert(sscanf(cval.c_str(),"%lf %lf",&r,&i)==2);
    f = std::complex<double>(r,i);
  }

  void parse(std::complex<float>& f, const std::string& cval) {
    float r,i;
    assert(sscanf(cval.c_str(),"%f %f",&r,&i)==2);
    f = std::complex<float>(r,i);
  }

  template<class T>
    void get(const char* name, std::vector<T>& v) {
    int i = 0;
    v.resize(0);
    while (true) {
      char buf[4096];
      sprintf(buf,"%s[%d]",name,i++);
      if (!has(buf))
	break;
      T val;
      parse(val,get(buf));
      std::cout << Grid::GridLogMessage << loghead() << "Set " << buf << " to " << val << std::endl;
      v.push_back(val);
    }
  }

  template<class T>
    void get(const char* name, T& f) {
    parse(f,get(name));
    std::cout << Grid::GridLogMessage << loghead() << "Set " << name << " to " << f << std::endl;
  }

  
};
