#include "characters.h"
#include "slint.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <SDL2/SDL.h>
#include <cmath>

using namespace std;

int parity = 0;
string parity_string("even");
string datapath("/files/math/projects/character-sum-data/max/");

typedef struct {
    unsigned int number;
    unsigned int max_location;
    double absmax;
    double argmax;
    double arggauss;
} chardata;


void rational_in_interval(double start, double end, int &x, int &y) {
    int a = 0;
    int b = 1;
    int c = 1;
    int d = 1;

    if(start == 0) {
        if(1.0/end < ceil(1/end)) {
            x = 1;
            y = (int)(ceil(1/end));
            return;
        }
        else {
            x = 1;
            y = (int)(ceil(1/end) + 1);
        }
    }

    double z = (a + c)/(double)(b + d);
    while(1) {
        if(start < z && z < end) {
            x = a+c;
            y = b+d;
            return;
        }
        else if(start < z) {
            c = a+c;
            d = b+d;
        }
        else {
            a = a+c;
            b = b+d;
        }
        z = (a+c)/(double)(b+d);
    }
}

void chi_with_denominator(int q, vector<int> &nlist, int D) {
    string padding;
    int qq = q;
    if(qq < 10)         padding = "0000000";
    else if(qq < 100)    padding = "000000";
    else if(qq < 1000)    padding = "00000";
    else if(qq < 10000)    padding = "0000";
    else if(qq < 100000)    padding = "000";
    else if(qq < 1000000)    padding = "00";
    else if(qq < 10000000)    padding = "0";
    else padding = "";

    string infilename = datapath + "/" + parity_string + "/" + padding + to_string(q);
    ifstream infile(infilename);

    chardata data;

    infile.read( (char*)&data, sizeof(data));
    int a,b;
    while(!infile.eof()) {
        rational_in_interval(data.max_location/(double)q, (data.max_location + 1)/(double)q, a, b);
        if(b == D) {
            nlist.push_back(data.number);
        }
        infile.read( (char*)&data, sizeof(data));
    }
}



class GraphingSurface {
public:
    double left;
    double right;
    double top;
    double bottom;

    int height;
    int width;

    SDL_Window * window = nullptr;
    SDL_Renderer * renderer = nullptr;

    SDL_Point * last_lines;
    SDL_Point * next_lines;
    int count = 0;
    bool first_update = true;
    int delay = 1000;
    
    GraphingSurface(int _width,
                    int _height,
                    double _left,
                    double _right,
                    double _bottom,
                    double _top) {

        height = _height;
        width = _width;
        left = _left;
        right = _right;
        top = _top;
        bottom = _bottom;

        last_lines = new SDL_Point[delay];
        next_lines = new SDL_Point[delay];
        next_lines[0].x = convertx(0);
        next_lines[0].y = converty(0);
        count = 1;


        window = SDL_CreateWindow("Hello World!", 100, 100, height, width, SDL_WINDOW_SHOWN | SDL_WINDOW_FULLSCREEN_DESKTOP);
        if (window == nullptr){
            std::cout << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
            SDL_Quit();
            exit(1);
        }

        //renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
        if (renderer == nullptr) {
            SDL_DestroyWindow(window);
            std::cout << "SDL_CreateRenderer Error: " << SDL_GetError() << std::endl;
            SDL_Quit();
            exit(1);
        }

    }

    int convertx(double x) {
        return (x - left) * width/(right - left);
    }
    int converty(double y) {
        return (top - y) * height/(top - bottom);
    }

    void line(double x1, double y1, double x2, double y2) {
        int X1 = (x1 - left) * width/(right - left);
        int X2 = (x2 - left) * width/(right - left);
        int Y1 = (top - y1) * height/(top - bottom);
        int Y2 = (top - y2) * height/(top - bottom);

        next_lines[count].x = X2;
        next_lines[count].y = Y2;
        count++;
        if(count == delay) update();
        //SDL_RenderDrawLine(renderer, X1, Y1, X2, Y2);
    }
    void update() {
        if(!first_update) {
            SDL_SetRenderDrawColor(renderer, 255,255,255,255);
            SDL_RenderDrawLines(renderer, last_lines, delay);
        }
        first_update = false;
        SDL_SetRenderDrawColor(renderer, 255,0,0,255);
        SDL_RenderDrawLines(renderer, next_lines, delay);
        SDL_RenderPresent(renderer);
        SDL_Point * tmp = next_lines;
        next_lines = last_lines;
        last_lines = tmp;
        count = 0;
        SDL_Delay(15);
    }
    void clear() {
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        first_update = true;
        next_lines[0].x = convertx(0);
        next_lines[0].y = converty(0);
        count = 1;
    }

    void change_color() {
        SDL_SetRenderDrawColor(renderer, rand() % 255, rand() % 255, rand() % 255, 255);
    }

    ~GraphingSurface() {
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
    }

    //http://stackoverflow.com/questions/20233469/how-do-i-take-and-save-a-bmp-screenshot-in-sdl-2
    bool saveBMP(std::string filepath) {
        SDL_Surface* saveSurface = NULL;
        SDL_Surface* infoSurface = NULL;
        infoSurface = SDL_GetWindowSurface(window);
        if (infoSurface == NULL) {
            std::cerr << "Failed to create info surface from window in saveScreenshotBMP(string), SDL_GetError() - " << SDL_GetError() << "\n";
        } else {
            unsigned char * pixels = new (std::nothrow) unsigned char[infoSurface->w * infoSurface->h * infoSurface->format->BytesPerPixel];
            if (pixels == 0) {
                std::cerr << "Unable to allocate memory for screenshot pixel data buffer!\n";
                return false;
            } else {
                if (SDL_RenderReadPixels(renderer, &infoSurface->clip_rect, infoSurface->format->format, pixels, infoSurface->w * infoSurface->format->BytesPerPixel) != 0) {
                    std::cerr << "Failed to read pixel data from SDL_Renderer object. SDL_GetError() - " << SDL_GetError() << "\n";
                    pixels = NULL;
                    return false;
                } else {
                    saveSurface = SDL_CreateRGBSurfaceFrom(pixels, infoSurface->w, infoSurface->h, infoSurface->format->BitsPerPixel, infoSurface->w * infoSurface->format->BytesPerPixel, infoSurface->format->Rmask, infoSurface->format->Gmask, infoSurface->format->Bmask, infoSurface->format->Amask);
                    if (saveSurface == NULL) {
                        std::cerr << "Couldn't create SDL_Surface from renderer pixel data. SDL_GetError() - " << SDL_GetError() << "\n";
                        return false;
                    }
                    SDL_SaveBMP(saveSurface, filepath.c_str());
                    SDL_FreeSurface(saveSurface);
                    saveSurface = NULL;
                }
                delete[] pixels;
            }
            SDL_FreeSurface(infoSurface);
            infoSurface = NULL;
        }
        return true;
    }
};

void pause_animation() {
    SDL_Event e;
    bool paused = true;
    while(paused) {
        SDL_WaitEvent(&e);
        if(e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_SPACE)
            paused = false;
    }
}

int main(int argc, char ** argv) {
    int q = atoi(argv[1]);
    int D = 1;
    if(argc > 2) D = atoi(argv[2]);
    if(argc > 3) {
        parity_string = argv[3];
        if(parity_string == "odd") {
            parity = 1;
        }
    }
    
    srand(time(NULL));
    if (SDL_Init(SDL_INIT_EVERYTHING) != 0){
        std::cerr << "SDL_Init Error: " << SDL_GetError() << std::endl;
        exit(1);
    }
    SDL_DisplayMode mode;
    SDL_GetCurrentDisplayMode(0, &mode);
    cout << mode.w << " " << mode.h << endl;
    double ratio = mode.w/(double)mode.h;
    GraphingSurface surface(mode.w, mode.h, -1*ratio, 1*ratio, -1, 1);
    DirichletGroup G(q);
    vector<int> nlist;
    if(D > 1) {
        chi_with_denominator(q, nlist, D);
    }
    SDL_Event e;
    bool quit = false;

    int update_delay = 500;
    int RR = 255;
    int GG = 255;
    int BB = 255;
    int number_of_characters = q;
    if(D > 1) number_of_characters = nlist.size();
    for(int k = 2; k < number_of_characters && !quit; k++) {
        if(GCD(k,q) != 1) continue;
        DirichletCharacter chi;
        if(D > 1)
            chi = G.character(nlist[k]);
        else
            chi = G.character(k);
        if(!chi.is_primitive()) continue;
        complex<double> M = abs(chi.max(NULL)) * chi.gauss_sum()/sqrt(q);
        complex<double> S = 0.0;
        int end = q;
        if(chi.is_even() == parity) continue;
        if(parity == 1) end = q/2;
        //SDL_Delay(1000);
        for(int n = 0; n < end && !quit; n++) {
            complex<double> S2 = S + chi.value(n)/M;
            surface.line(S.real(), S.imag(), S2.real(), S2.imag());
            S = S2;
            //cout << n << endl;
            while (SDL_PollEvent(&e)){
		if(e.type == SDL_QUIT){
                    quit = true;
		}
		if(e.type == SDL_KEYDOWN) {
                    switch(e.key.keysym.sym) {
                        case SDLK_RIGHT:
                            n = end;
                            break;
                        case SDLK_LEFT:
                            n = end;
                            k -= 3;
                            break;
                        case SDLK_p:
                            update_delay += 10;
                            break;
                        case SDLK_o:
                            update_delay -= 10;
                            update_delay = max(update_delay, 1);
                            break;
                        case SDLK_SPACE:
                            pause_animation();
                            break;
                        case SDLK_q:
                            quit = true;
                            break;
                        case SDLK_r:
                            RR = 255; GG = 0; BB = 0;
                            SDL_SetRenderDrawColor(surface.renderer, RR, GG, BB, 255);
                            break;
                        case SDLK_g:
                            RR = 0; GG = 255; BB = 0;
                            SDL_SetRenderDrawColor(surface.renderer, RR, GG, BB, 255);
                            break;
                        case SDLK_b:
                            RR = 0; GG = 0; BB = 255;
                            SDL_SetRenderDrawColor(surface.renderer, RR, GG, BB, 255);
                            break;
                        case SDLK_s:
                            surface.saveBMP("image.bmp");
                            break;
                    }
                }
            }
            //if(n % update_delay == 0) surface.update();
        }
        surface.update();
        SDL_Delay(200);
        surface.clear();
        //SDL_SetRenderDrawColor(surface.renderer, 0, 0, 0, 255);
        //SDL_RenderClear(surface.renderer);
        //SDL_SetRenderDrawColor(surface.renderer, RR, GG, BB, 255);
        //surface.change_color();
    }

    //SDL_Delay(2000);

    return 0;
}
