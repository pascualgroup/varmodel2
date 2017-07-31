#ifndef __varmodel_hpp__
#define __varmodel_hpp__

namespace varmodel {

void initialize();

void load_checkpoint();
void save_checkpoint();

void run();

void update_biting_rate();

}; // namespace varmodel

#endif /* __varmodel_hpp__ */
