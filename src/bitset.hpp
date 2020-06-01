#pragma once

#include <cstdint>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <random>
#include <limits>
#include <cassert>

#define BITS 64

namespace triangulator {
// FIXEDBITSET
template <size_t chunks>
class FBitset {
 public:
  uint64_t data_[chunks];
  void Clear() {
    for (size_t i=0;i<chunks;i++){
      data_[i] = 0;
    }
  }
  FBitset() {
    Clear();
  }
  bool operator<(const FBitset<chunks>& other) const {
    for (size_t i=0;i<chunks;i++){
      if (data_[i]<other.data_[i]) return true;
      else if(data_[i]>other.data_[i]) return false;
    }
    return false;
  }
  bool operator==(const FBitset<chunks>& other) const {
    for (size_t i=0;i<chunks;i++){
      if (data_[i] != other.data_[i]) return false;
    }
    return true;
  }
  bool operator!=(const FBitset<chunks>& other) const {
    for (size_t i=0;i<chunks;i++){
      if (data_[i] != other.data_[i]) return true;
    }
    return false;
  }
  FBitset<chunks> operator|(const FBitset<chunks>& other) const {
    FBitset<chunks> ret;
    for (size_t i=0;i<chunks;i++){
      ret.data_[i] = data_[i] | other.data_[i];
    }
    return ret;
  }
  FBitset<chunks> operator&(const FBitset<chunks>& other) const {
    FBitset<chunks> ret;
    for (size_t i=0;i<chunks;i++){
      ret.data_[i] = data_[i] & other.data_[i];
    }
    return ret;
  }
  FBitset<chunks> operator~() const {
    FBitset<chunks> ret;
    for (size_t i=0;i<chunks;i++){
      ret.data_[i] = (~data_[i]);
    }
    return ret;
  }
  void Set(size_t i, bool v) {
    if (v) {
      data_[i/BITS] |= ((uint64_t)1 << (uint64_t)(i%BITS));
    } else {
      data_[i/BITS] &= (~((uint64_t)1 << (uint64_t)(i%BITS)));
    }
  }
  void SetTrue(size_t i) {
    data_[i/BITS] |= ((uint64_t)1 << (uint64_t)(i%BITS));
  }
  void SetFalse(size_t i) {
    data_[i/BITS] &= (~((uint64_t)1 << (uint64_t)(i%BITS)));
  }
  void SetTrue(const std::vector<size_t>& v) {
    for (size_t x : v) {
      SetTrue(x);
    }
  }
  void SetTrue(const std::vector<int>& v) {
    for (int x : v) {
      SetTrue(x);
    }
  }
  void SetFalse(const std::vector<int>& v) {
    for (int x : v) {
      SetFalse(x);
    }
  }
  void FillTrue() {
    for (size_t i=0;i<chunks;i++){
      data_[i] = ~0;
    }
  }
  void FillUpTo(size_t n) {
    for (size_t i=0;i<chunks;i++){
      if ((i+1)*BITS <= n) {
        data_[i] = ~0;
      } else if (i*BITS < n) {
        for (size_t j=i*BITS;j<n;j++){
          SetTrue(j);
        }
      } else {
        return;
      }
    }
  }
  bool Get(size_t i) const {
    return data_[i/BITS] & ((uint64_t)1 << (uint64_t)(i%BITS));
  }
  bool IsEmpty() const {
    for (size_t i=0;i<chunks;i++){
      if (data_[i]) return false;
    }
    return true;
  }
  void operator |= (const FBitset<chunks>& rhs) {
    for (size_t i=0;i<chunks;i++){
      data_[i] |= rhs.data_[i];
    }
  }
  void operator &= (const FBitset<chunks>& rhs) {
    for (size_t i=0;i<chunks;i++){
      data_[i] &= rhs.data_[i];
    }
  }
  void TurnOff(const FBitset<chunks>& rhs) {
    for (size_t i=0;i<chunks;i++){
      data_[i] &= (~rhs.data_[i]);
    }
  }
  void InvertAnd(const FBitset<chunks>& rhs) {
    for (size_t i=0;i<chunks;i++){
      data_[i] = (~data_[i]) & rhs.data_[i];
    }
  }
  void SetNeg(const FBitset<chunks>& rhs) {
    for (size_t i=0;i<chunks;i++){
      data_[i] = ~rhs.data_[i];
    }
  }
  void SetNegAnd(const FBitset<chunks>& rhs1, const FBitset<chunks>& rhs2) {
    for (size_t i=0;i<chunks;i++){
      data_[i] = (~rhs1.data_[i]) & rhs2.data_[i];
    }
  }
  void SetAnd(const FBitset<chunks>& rhs1, const FBitset<chunks>& rhs2) {
    for (size_t i=0;i<chunks;i++){
      data_[i] = rhs1.data_[i] & rhs2.data_[i];
    }
  }
  bool Subsumes(const FBitset<chunks>& other) const {
    for (size_t i=0;i<chunks;i++){
      if ((data_[i] | other.data_[i]) != data_[i]) return false;
    }
    return true;
  }
  std::vector<int> Elements() const {
    std::vector<int> ret;
    for (size_t i=0;i<chunks;i++){
      uint64_t td = data_[i];
      while (td) {
        ret.push_back(i*BITS + __builtin_ctzll(td));
        td &= ~-td;
      }
    }
    return ret;
  }
  int Popcount() const {
    int cnt = 0;
    for (size_t i=0;i<chunks;i++) {
      cnt += __builtin_popcountll(data_[i]);
    }
    return cnt;
  }
  bool Intersects(const FBitset<chunks>& other) const {
    for (size_t i=0;i<chunks;i++){
      if (data_[i] & other.data_[i]) return true;
    }
    return false;
  }
  int IntersectionPopcount(const FBitset<chunks>& other) const {
    int cnt = 0;
    for (size_t i=0;i<chunks;i++) {
      cnt += __builtin_popcountll(data_[i] & other.data_[i]);
    }
    return cnt;
  }
  int First() const {
    for (size_t i=0;i<chunks;i++) {
      if (data_[i]) {
        return i*BITS + __builtin_ctzll(data_[i]);
      }
    }
    return chunks * BITS;
  }

  class FBitsetIterator {
   private:
    const FBitset<chunks>* const bitset_;
    size_t pos_;
    uint64_t tb_;
   public:
    FBitsetIterator(const FBitset<chunks>* const bitset, size_t pos, uint64_t tb) : bitset_(bitset), pos_(pos), tb_(tb) { }
    bool operator!=(const FBitsetIterator& other) const {
      return pos_ != other.pos_ || tb_ != other.tb_;
    }
    const FBitsetIterator& operator++() {
      tb_ &= ~-tb_;
      while (tb_ == 0 && pos_ < chunks) {
        pos_++;
        if (pos_ < chunks) {
          tb_ = bitset_->data_[pos_];
        }
      }
      return *this;
    }
    int operator*() const {
      return pos_*BITS + __builtin_ctzll(tb_);
    }
  };

  FBitsetIterator begin() const {
    size_t pos = 0;
    while (pos < chunks && data_[pos] == 0) {
      pos++;
    }
    if (pos < chunks) {
      return FBitsetIterator(this, pos, data_[pos]);
    } else {
      return FBitsetIterator(this, pos, 0);
    }
  }
  FBitsetIterator end() const {
    return FBitsetIterator(this, chunks, 0);
  }
};

template <size_t chunks>
class FBitsetSet {
 public:
  FBitsetSet() {}
  FBitsetSet(size_t capacity, double load_factor) {
    load_factor_ = load_factor;
    assert(chunks > 0);
    assert(load_factor_ >= 1.1);
    capacity_ = NextPrime((capacity + 1) * load_factor_);
    assert((size_t)(capacity_ * load_factor_) > capacity_);
    container_.resize(capacity_ * chunks);
  }
  bool Insert(const FBitset<chunks>& bitset) {
    size_t ind = Hash(bitset.data_, capacity_);
    while (1) {
      if (Zero(IndToPtr(ind, container_))) break;
      else if (Equal(IndToPtr(ind, container_), bitset.data_)) return false;
      else {
        ind++;
        if (ind == capacity_) {
          ind = 0;
        }
      }
    }
    Copy(bitset.data_, IndToPtr(ind, container_));
    elements_++;
    if ((size_t)(elements_ * load_factor_) > capacity_) {
      Resize();
      assert((size_t)(elements_ * load_factor_) < capacity_);
    }
    return true;
  }
  bool Contains(const FBitset<chunks>& bitset) const {
    size_t ind = Hash(bitset.data_, capacity_);
    while (1) {
      if (Zero(IndToPtr(ind, container_))) return false;
      else if (Equal(IndToPtr(ind, container_), bitset.data_)) return true;
      else {
        ind++;
        if (ind == capacity_) {
          ind = 0;
        }
      }
    }
  }
  bool Inited() const {
    return capacity_ > 0;
  }
  size_t ContainerSize() const {
    return container_.size();
  }
  std::vector<FBitset<chunks>> Vector() const {
    std::vector<FBitset<chunks>> ret;
    ret.reserve(elements_);
    for (size_t i = 0; i < capacity_; i++) {
      if (Zero(IndToPtr(i, container_))) continue;
      ret.push_back(FBitset<chunks>());
      for (size_t j = 0; j < chunks; j++) {
        ret.back().data_[j] = container_[i*chunks + j];
      }
    }
    assert(ret.size() == elements_);
    return ret;
  }
 private:
  size_t elements_ = 0;
  double load_factor_ = 0;
  size_t capacity_ = 0;
  std::vector<uint64_t> container_;

  bool IsPrime(size_t n) const {
    if (n < 2) return false;
    if (n < 4) return true;
    if (n%2 == 0) return false;
    for (size_t i=3;i*i<=n;i+=2) {
      if (n%i == 0) return false;
    }
    return true;
  }

  size_t NextPrime(size_t n) const {
    while (!IsPrime(n)) {
      n++;
    }
    return n;
  }

  size_t Hash(const uint64_t* const data, size_t mod) const {
    uint64_t h = 0;
    for (size_t i=0;i<chunks;i++){
      h = h*65599 + data[i];
    }
    h %= (uint64_t)mod;
    return h;
  }

  bool Equal(const uint64_t* const data1, const uint64_t* const data2) const {
    for (size_t i=0;i<chunks;i++){
      if (data1[i] != data2[i]) return false;
    }
    return true;
  }

  bool Zero(const uint64_t* const data) const {
    for (size_t i=0;i<chunks;i++){
      if (data[i]) return false;
    }
    return true;
  }

  const uint64_t* IndToPtr(size_t ind, const std::vector<uint64_t>& vec) const {
    return vec.data() + (ind * chunks);
  }

  uint64_t* IndToPtr(size_t ind, std::vector<uint64_t>& vec) const {
    return vec.data() + (ind * chunks);
  }

  void Copy(const uint64_t* const from, uint64_t* const to) {
    for (size_t i=0;i<chunks;i++){
      to[i] = from[i];
    }
  }

  void Resize() {
    size_t new_capacity = NextPrime((capacity_ * 3) / 2);
    std::vector<uint64_t> new_container(new_capacity * chunks);
    for (size_t i = 0; i < capacity_; i++) {
      if (Zero(IndToPtr(i, container_))) continue;
      size_t ind = Hash(IndToPtr(i, container_), new_capacity);
      while (!Zero(IndToPtr(ind, new_container))) {
        ind++;
        if (ind == new_capacity) {
          ind = 0;
        }
      }
      Copy(IndToPtr(i, container_), IndToPtr(ind, new_container));
    }
    container_ = std::move(new_container);
    capacity_ = new_capacity;
  }
};

template<size_t chunks>
class FBitsetMap {
 public:
  FBitsetMap() {}
  FBitsetMap(size_t capacity, double load_factor) {
    load_factor_ = load_factor;
    assert(chunks > 0);
    assert(load_factor_ >= 1.1);
    capacity_ = NextPrime((capacity + 1) * load_factor_);
    assert((size_t)(capacity_ * load_factor_) > capacity_);
    container_.resize(capacity_ * chunks);
    values_.resize(capacity_);
  }
  std::pair<int, bool> Insert(const FBitset<chunks>& bitset, int value, bool replace=false) {
    assert(value > 0);
    size_t ind = Hash(bitset.data_, capacity_);
    while (1) {
      if (Zero(IndToPtr(ind, container_))) break;
      else if (Equal(IndToPtr(ind, container_), bitset.data_)) {
        assert(values_[ind] > 0);
        if (replace) {
          values_[ind] = value;
        }
        return {values_[ind], false};
      } else {
        ind++;
        if (ind == capacity_) {
          ind = 0;
        }
      }
    }
    Copy(bitset.data_, IndToPtr(ind, container_));
    elements_++;
    values_[ind] = value;
    if ((size_t)(elements_ * load_factor_) > capacity_) {
      Resize();
      assert((size_t)(elements_ * load_factor_) < capacity_);
    }
    return {value, true};
  }
  int Get(const FBitset<chunks>& bitset, bool expect) const {
    size_t ind = Hash(bitset.data_, capacity_);
    while (1) {
      if (Zero(IndToPtr(ind, container_))) {
        assert(!expect);
        return 0;
      } else if (Equal(IndToPtr(ind, container_), bitset.data_)) {
        assert(values_[ind] > 0);
        return values_[ind];
      } else {
        ind++;
        if (ind == capacity_) {
          ind = 0;
        }
      }
    }
    assert(0);
  }
  bool Inited() const {
    return capacity_ > 0;
  }
  size_t ContainerSize() const {
    return container_.size();
  }
 private:
  size_t elements_ = 0;
  double load_factor_ = 0;
  size_t capacity_ = 0;
  std::vector<uint64_t> container_;
  std::vector<int> values_;

  bool IsPrime(size_t n) const {
    if (n < 2) return false;
    if (n < 4) return true;
    if (n%2 == 0) return false;
    for (size_t i=3;i*i<=n;i+=2) {
      if (n%i == 0) return false;
    }
    return true;
  }

  size_t NextPrime(size_t n) const {
    while (!IsPrime(n)) {
      n++;
    }
    return n;
  }

  size_t Hash(const uint64_t* const data, size_t mod) const {
    uint64_t h = 0;
    for (size_t i=0;i<chunks;i++){
      h = h*65599 + data[i];
    }
    h %= (uint64_t)mod;
    return h;
  }

  bool Equal(const uint64_t* const data1, const uint64_t* const data2) const {
    for (size_t i=0;i<chunks;i++){
      if (data1[i] != data2[i]) return false;
    }
    return true;
  }

  bool Zero(const uint64_t* const data) const {
    for (size_t i=0;i<chunks;i++){
      if (data[i]) return false;
    }
    return true;
  }

  const uint64_t* IndToPtr(size_t ind, const std::vector<uint64_t>& vec) const {
    return vec.data() + (ind * chunks);
  }

  uint64_t* IndToPtr(size_t ind, std::vector<uint64_t>& vec) const {
    return vec.data() + (ind * chunks);
  }

  void Copy(const uint64_t* const from, uint64_t* const to) {
    for (size_t i=0;i<chunks;i++){
      to[i] = from[i];
    }
  }

  void Resize() {
    size_t new_capacity = NextPrime((capacity_ * 3) / 2);
    std::vector<uint64_t> new_container(new_capacity * chunks);
    std::vector<int> new_values(new_capacity);
    for (size_t i = 0; i < capacity_; i++) {
      if (Zero(IndToPtr(i, container_))) continue;
      size_t ind = Hash(IndToPtr(i, container_), new_capacity);
      while (!Zero(IndToPtr(ind, new_container))) {
        ind++;
        if (ind == new_capacity) {
          ind = 0;
        }
      }
      Copy(IndToPtr(i, container_), IndToPtr(ind, new_container));
      assert(values_[i] > 0);
      new_values[ind] = values_[i];
    }
    container_ = std::move(new_container);
    values_ = std::move(new_values);
    capacity_ = new_capacity;
  }
};


template<size_t chunks>
class FLBSieve {
 public:
  FLBSieve() {}
  FLBSieve(size_t len, int k) : len_(len) {
    assert(len_ <= chunks * BITS);
    assert(k>=0);
    containers_.resize(k+1);
    masks_.resize(k+1);
    elements_.resize(k+1);
    for (int i=0;i<=k;i++){
      containers_[i].push_back({});
    }
  }
  void Insert(const FBitset<chunks>& bs, int lb) {
    assert(lb < masks_.size() && lb>=0);
    maxlb_ = std::max(maxlb_, lb);
    int mask = GetMask(lb, bs);
    for (size_t i=0;i<chunks;i++) {
      containers_[lb][mask].push_back(bs.data_[i]);
    }
    elements_[lb]++;
    if (elements_[lb] > (1 << (2*masks_[lb].size()))) {
      Resize(lb);
      assert(elements_[lb] <= (1 << (2*masks_[lb].size())));
    }
  }
  int Get(const FBitset<chunks>& bs, int k) {
    for (int i=k;i<=maxlb_;i++) {
      int mask = GetMask(i, bs);
      if (GetCont(i, bs, 0)) {
        return i;
      }
      for (int sub=0;(sub=(sub-mask)&mask);){
        if (GetCont(i, bs, sub)) {
          return i;
        }
      }
    }
    return 0;
  }
  size_t TotElements() const {
    size_t ret = 0;
    for (size_t e : elements_) {
      ret += e;
    }
    return ret;
  }
 private:
  int maxlb_ = 0;
  size_t len_ = 0;
  std::vector<size_t> elements_;
  std::vector<std::vector<FBitset<chunks>>> masks_;
  std::vector<std::vector<std::vector<uint64_t>>> containers_;
  int GetMask(int lb, const FBitset<chunks>& bs) {
    int mask = 0;
    for (int i=0;i<(int)masks_[lb].size();i++) {
      if (bs.Intersects(masks_[lb][i])) {
        mask |= (1<<i);
      }
    }
    return mask;
  }
  bool GetCont(int lb, const FBitset<chunks>& bs, int cont) {
    size_t cs = containers_[lb][cont].size()/chunks;
    for (size_t i=0;i<cs;i++) {
      bool issub = true;
      for (size_t j=0;j<chunks;j++) {
        if ((bs.data_[j] | containers_[lb][cont][i*chunks + j]) != bs.data_[j]) {
          issub = false;
          break;
        }
      }
      if (issub) {
        return true;
      }
    }
    return false;
  }
  int NIntersect(int lb, const FBitset<chunks>& bs) {
    int n=0;
    for (int cont=0;cont<(int)containers_[lb].size();cont++) {
      size_t cs = containers_[lb][cont].size()/chunks;
      for (size_t i=0;i<cs;i++) {
        bool iss = false;
        for (size_t j=0;j<chunks;j++) {
          if (bs.data_[j] & containers_[lb][cont][i*chunks + j]) {
            iss = true;
            break;
          }
        }
        if (iss) {
          n++;
        }
      }
    }
    return n;
  }
  void Resize(int lb) {
    int nmasks = masks_[lb].size()+1;
    //std::cerr<<"resize "<<lb<<" "<<nmasks<<" "<<elements_[lb]<<std::endl;
    masks_[lb].clear();
    int ib = len_/nmasks;
    assert(nmasks*ib <= len_);
    assert(ib >= 2);
    int sp = 0;
    FBitset<chunks> mask;
    for (int i=0;i<nmasks;i++) {
      int ep = (i+1)*ib;
      assert(sp+ib <= ep);
      assert(ep<=len_);
      mask.Clear();
      for (int j=sp;j<ep;j++){
        mask.SetTrue(j);
        int ic = NIntersect(lb, mask);
        if (ic >= elements_[lb]/2 || j+1 == ep) {
          //std::cerr<<"nmask "<<mask.Popcount()<<" "<<ic<<std::endl;
          masks_[lb].push_back(mask);
          sp=ep;
          break;
        }
      }
    }
    assert((int)masks_[lb].size() == nmasks);
    std::vector<int> cnts(1 << masks_[lb].size());
    for (int cont=0;cont<(int)containers_[lb].size();cont++) {
      size_t cs = containers_[lb][cont].size()/chunks;
      for (int i=0;i<cs;i++) {
        int tmask = 0;
        for (int k=0;k<(int)masks_[lb].size();k++){
          bool iss = false;
          for (int j=0;j<chunks;j++) {
            if (masks_[lb][k].data_[j] & containers_[lb][cont][i*chunks + j]) {
              iss = true;
              break;
            }
          }
          if (iss) {
            tmask |= (1<<k);
          }
        }
        cnts[tmask]++;
      }
    }
    std::vector<std::vector<uint64_t>> new_containers(cnts.size());
    for (int i=0;i<(int)cnts.size();i++){
      new_containers[i].reserve(chunks*cnts[i]*2);
    }
    for (int cont=0;cont<(int)containers_[lb].size();cont++) {
      size_t cs = containers_[lb][cont].size()/chunks;
      for (int i=0;i<cs;i++) {
        int tmask = 0;
        for (int k=0;k<(int)masks_[lb].size();k++){
          bool iss = false;
          for (int j=0;j<chunks;j++) {
            if (masks_[lb][k].data_[j] & containers_[lb][cont][i*chunks + j]) {
              iss = true;
              break;
            }
          }
          if (iss) {
            tmask |= (1<<k);
          }
        }
        for (int j=0;j<chunks;j++) {
          new_containers[tmask].push_back(containers_[lb][cont][i*chunks + j]);
        }
      }
    }
    containers_[lb] = new_containers;
    for (int i=0;i<(int)cnts.size();i++) {
      assert(containers_[lb][i].size() == cnts[i]*chunks);
    }
  }
};

} // namespace triangulator