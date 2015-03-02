#pragma once

#include "MapShaper.h"

/**
 * Shapes the map into the Red Blob Games logo
 */
class BlobMapShaper: public MapShaper {
	public:
		BlobMapShaper();
		virtual bool isLand(glm::vec2);
};
