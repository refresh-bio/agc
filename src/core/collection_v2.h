#ifndef _COLLECTION_V2_H
#define _COLLECTION_V2_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include "collection.h"
#include "collection_v1.h"

class CCollection_V2 : public CCollection_V1
{

public:
	CCollection_V2() : CCollection_V1() {}
	virtual ~CCollection_V2() {};

	void serialize(vector<uint8_t>& data_main, vector<vector<uint8_t>>& data_details, bool store_date_time, uint32_t _details_batch_size);
	bool deserialize_main(vector<uint8_t>& data_main, bool create_maps);
	bool deserialize_details(vector<uint8_t>& zstd_data_details, size_t raw_size, bool deserialize_details);
};

// EOF
#endif