#include <LibIV/libiv.h>
#include <LibIV/CppHelpers.h>
#include <highgui.h>

namespace Poisson
{
	void possionImageEditing(TensorV3i & src_pixel, IplImage * dst_img);

	struct INFO
	{
		int N;    // The number of neighbors
		v3i self; // (x,y,index)
		v4i left; // (isExist, x, y, index)
		v4i right;
		v4i up;
		v4i down;
	};

	inline double PIE_b(const INFO & pix, const TensorInt & src, const TensorInt & dst)
	{
		double ans = 0;

		int width = src.cols();
		int height = src.rows();

		ans += pix.N * src(pix.self[0],pix.self[1]);

		/*if(pix.left[0])
		{*/
		ans -= src(pix.left[1],pix.left[2]);

		if(pix.left[1] == 0)
			ans += src(pix.self[0],pix.self[1]);
		else if(pix.left[3] == -1)
			ans += dst(pix.left[1],pix.left[2]);
		//}

		/*if(pix.right[0])
		{*/
		ans -= src(pix.right[1],pix.right[2]);

		if(pix.right[1] == width - 1)
			ans += src(pix.self[0],pix.self[1]);
		else if(pix.right[3] == -1)
			ans += dst(pix.right[1],pix.right[2]);
		//}

		/*if(pix.up[0])
		{*/
		ans -= src(pix.up[1],pix.up[2]);

		if(pix.up[2] == 0)
			ans += src(pix.self[0],pix.self[1]);
		else if(pix.up[3] == -1)
			ans += dst(pix.up[1],pix.up[2]);
		//}

		/*if(pix.down[0])
		{*/
		ans -= src(pix.down[1],pix.down[2]);

		if(pix.down[2] == height-1)
			ans += src(pix.self[0],pix.self[1]);
		else if(pix.down[3] == -1)
			ans += dst(pix.down[1],pix.down[2]);
		//}

		return ans;
	}
}

