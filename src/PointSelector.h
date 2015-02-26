#pragma once

#include <vector>
#include <glm/glm.hpp>

/**
 * Abstract class used to select points to generate a map
 */
class PointSelector {
	protected:
		// The width of the map to generate
		int width;
		// The height of the map to generate
		int height;
	public:
		/**
		 * @param in int mapWidth The width of the map to generate
		 * @param in int mapHeight The height of the map to generate
		 */
		PointSelector(int mapWidth, int mapHeight){
			width = mapWidth;
			height = mapHeight;
		}

		/**
		 * Selects `number` of points for this map
		 *
		 * @param in int number The number of points to return
		 * @return Vector<glm::vec2> The selected points
		 */
		virtual std::vector<glm::vec2> select(int number) = 0;
};
