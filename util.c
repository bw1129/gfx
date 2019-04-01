#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <io.h>
#include <direct.h>
#include <wtypes.h>
#include <wingdi.h> 
#include <time.h>

#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>

#include <EL2/eyelink.h>
#include <EL2/w32_exptsppt2.h>
#include <EL2/w32_demo.h>

#include "exp.h"
#include "maths.h"

extern int 
	Break_it;

#define CAL_COL			5
#define CAL_PTS			(CAL_COL)
#define CAL_REP			3

int UTIL_time_delta(int reset)
{
	static long old_time, this_time;
	int delta;
	this_time = SDL_GetTicks();
	delta = this_time - old_time;
	if(reset) old_time = this_time;
	return delta;
}
void UTIL_print_time(FILE *fp)
{
    struct tm *newtime;
	char am_pm[] = "AM";
	time_t long_time;

	time( &long_time );                /* Get time as long integer. */
	newtime = localtime( &long_time ); /* Convert to local time. */
	if( newtime->tm_hour > 12 )        /* Set up extension. */
		strcpy( am_pm, "PM" );
	if( newtime->tm_hour > 12 )        /* Convert from 24-hour */
	newtime->tm_hour -= 12;    /*   to 12-hour clock.  */
	if( newtime->tm_hour == 0 )        /*Set hour to 12 if midnight. */
	newtime->tm_hour = 12;
	printf( "%.19s %s\n", asctime( newtime ), am_pm );
	fprintf( fp, "%.19s %s\n", asctime( newtime ), am_pm );
}

void IO_toggle_grab(void)
{
	SDL_GrabMode mode;

	printf("Ctrl-G: toggling input grab!\n");
	mode = SDL_WM_GrabInput(SDL_GRAB_QUERY);
	if ( mode == SDL_GRAB_ON ) {
		printf("Grab was on\n");
	} else {
		printf("Grab was off\n");
	}
	mode = SDL_WM_GrabInput(!mode);
	if ( mode == SDL_GRAB_ON ) {
		printf("Grab is now on\n");
	} else {
		printf("Grab is now off\n");
	}
}

int IO_handle_event(SDL_Event *event, Uint8 type, SDLKey key, Uint16 *x, Uint16 *y, Uint8 *b)
{
	int done=0;
	if(event->type==SDL_KEYDOWN) {
		if(event->key.keysym.sym == SDLK_ESCAPE ) { 
					leave("Escape pressed.\n");
		}
	}
	if(event->type == SDL_KEYDOWN) {
		if(event->key.keysym.sym == SDLK_x) Break_it =1;
	}
	if(event->type == type) {
		switch(event->type ) {
			case SDL_KEYDOWN:
				if(event->key.keysym.sym == key) {
					done = 1;
				}
				break;
			case SDL_MOUSEBUTTONDOWN:
				*x = event->button.x;
				*y = event->button.y;
				*b = event->button.button;
				done = 1;
				break;
			case SDL_QUIT:
				leave("SDL_QUIT\n");
				break;
			default:
				break;
		}
	}
	return(done);
}
void IO_check_escape(void)
{
	SDL_Event event;
	SDL_PollEvent(&event);
	IO_handle_event(&event, SDL_KEYDOWN, SDLK_ESCAPE, NULL, NULL, NULL);
}
int IO_check_key(SDLKey key)
{
	SDL_Event event;
	SDL_PollEvent(&event);
	return (IO_handle_event(&event, SDL_KEYDOWN, key, NULL, NULL, NULL));
}
void IO_flush_events(void)
{
	SDL_Event event;
	while(SDL_PollEvent(&event));
}
void IO_wait_for_event(Uint8 type, SDLKey key, Uint16 *x, Uint16 *y, Uint8 *b)
{
	SDL_Event event;
	do {
		SDL_PollEvent(&event);
	} while (! IO_handle_event(&event, type, key, x, y, b));
}
SDLKey IO_get_key(void)
{
	SDL_Event event;
	do {
		SDL_PollEvent(&event);
	} while (event.type!=SDL_KEYDOWN);
	if(event.key.keysym.sym == SDLK_ESCAPE ) { 
		leave("Escape pressed.\n");
	}
	return event.key.keysym.sym;
}
SDLKey IO_ask_key(void)
{
	SDL_Event event;
	SDL_PollEvent(&event);
	if(event.type==SDL_KEYDOWN) {
		if(event.key.keysym.sym == SDLK_ESCAPE ) { 
			leave("Escape pressed.\n");
		}
		return event.key.keysym.sym;
	} else {
		return 0;
	}
}
SDLKey IO_ask_key_release(void)
{
	SDL_Event event;
	SDL_PollEvent(&event);
	if(event.type==SDL_KEYUP) {
		if(event.key.keysym.sym == SDLK_ESCAPE ) { 
			leave("Escape pressed.\n");
		}
		return event.key.keysym.sym;
	} else {
		return 0;
	}
}