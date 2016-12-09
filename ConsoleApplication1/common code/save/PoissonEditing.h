#pragma  once
#include <LibIV/LibIV.h>

namespace PoissonEditing
{
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

		ans += pix.N * src(pix.self[0],pix.self[1]);


		if(pix.left[0])
		{
			ans -= src(pix.left[1],pix.left[2]);

			if(pix.left[3] == -1)
				ans += dst(pix.left[1],pix.left[2]);
		}

		if(pix.right[0])
		{
			ans -= src(pix.right[1],pix.right[2]);

			if(pix.right[3] == -1)
				ans += dst(pix.right[1],pix.right[2]);
		}

		if(pix.up[0])
		{
			ans -= src(pix.up[1],pix.up[2]);

			if(pix.up[3] == -1)
				ans += dst(pix.up[1],pix.up[2]);
		}

		if(pix.down[0])
		{
			ans -= src(pix.down[1],pix.down[2]);

			if(pix.down[3] == -1)
				ans += dst(pix.down[1],pix.down[2]);
		}

		return ans;
	}

	void poissonImageEditing_one_channel(TensorByte & ans, const TensorInt & src, const TensorInt & dst, 
		const TensorByte & msk);

	// =========================================================================
	// Poisson Image Editing
	// src_pixel: store the pixels to be composited into dst_img. All others are (-1,-1,-1).
	// dst_img:   the destination.

	void poissonImageEditing(TensorV3i & src_pixel, IplImage * dst_img);
}

