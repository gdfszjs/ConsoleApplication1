#pragma once

#include <string>
using namespace std;


namespace NMessage
{
	class NMessage
	{
	private:
		string msg;
	public:
		NMessage();
		NMessage(const string & m);

		const string & What();
	};


	NMessage::NMessage()
	{
		msg="Something wrong!";
	}

	NMessage::NMessage(const string & m) : msg(m) {};

	const string & NMessage::What()
	{
		return msg;
	}
}

