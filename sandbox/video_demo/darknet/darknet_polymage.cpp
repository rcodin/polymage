#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include "../simple_pool_allocator.h"
#include "../darknet.h"

#define isl_min(x,y) ((x) < (y) ? (x) : (y))
#define isl_max(x,y) ((x) > (y) ? (x) : (y))
#define isl_floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
extern "C" void  pipeline_opt(int  C, int  R, network *net, char* datacfg)
{
   // list *options = read_data_cfg(datacfg);
   //  char *name_list = option_find_str(options, "names", "data/names.list");
   //  char **names = get_labels(name_list);

   //  char *filename = "../messigray.png";
   //  float thresh = .5;
   //  float hier_thresh = .5;
   //  char *outfile = "try.png"; 
   //  int fullscreen = 0;


   //  image **alphabet = load_alphabet();
   //  // network *net = load_network(cfgfile, weightfile, 0);
   //  set_batch_network(net, 1);
   //  srand(2222222);
   //  double time;
   //  char buff[256];
   //  char *input = buff;
   //  int j;
   //  float nms=.3;
    
   //      if(filename){
   //          strncpy(input, filename, 256);
   //      } else {
   //          printf("Enter Image Path: ");
   //          fflush(stdout);
   //          input = fgets(input, 256, stdin);
   //          if(!input) return;
   //          strtok(input, "\n");
   //      }
   //      image im = load_image_color(input,0,0);
   //      image sized = letterbox_image(im, net->w, net->h);
   //      //image sized = resize_image(im, net->w, net->h);
   //      //image sized2 = resize_max(im, net->w);
   //      //image sized = crop_image(sized2, -((net->w - sized2.w)/2), -((net->h - sized2.h)/2), net->w, net->h);
   //      //resize_network(net, sized.w, sized.h);
   //      layer l = net->layers[net->n-1];

   //      box *boxes = (box*)calloc(l.w*l.h*l.n, sizeof(box));
   //      float **probs = (float**)calloc(l.w*l.h*l.n, sizeof(float *));
   //      for(j = 0; j < l.w*l.h*l.n; ++j) probs[j] = (float*)calloc(l.classes + 1, sizeof(float *));
   //      float **masks = 0;
   //      if (l.coords > 4){
   //          masks = (float**)calloc(l.w*l.h*l.n, sizeof(float*));
   //          for(j = 0; j < l.w*l.h*l.n; ++j) masks[j] = (float*)calloc(l.coords-4, sizeof(float *));
   //      }

   //      float *X = sized.data;
   //      time=what_time_is_it_now();
   //      network_predict(net, X);
   //      printf("%s: Predicted in %f seconds.\n", input, what_time_is_it_now()-time);
   //      get_region_boxes(l, im.w, im.h, net->w, net->h, thresh, probs, boxes, masks, 0, 0, hier_thresh, 1);
   //      //if (nms) do_nms_obj(boxes, probs, l.w*l.h*l.n, l.classes, nms);
   //      if (nms) do_nms_sort(boxes, probs, l.w*l.h*l.n, l.classes, nms);
   //      draw_detections(im, l.w*l.h*l.n, thresh, boxes, probs, masks, names, alphabet, l.classes);
   //      if(outfile){
   //          save_image(im, outfile);
   //      }
        

   //      free_image(im);
   //      free_image(sized);
   //      free(boxes);
   //      free_ptrs((void **)probs, l.w*l.h*l.n);
}


