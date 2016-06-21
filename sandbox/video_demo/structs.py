from cv2 import *
from common import clock, draw_str
from enum import Enum
import os.path
import ctypes

class AppType():
    NONE = -1
    UNSHARP = 0
    HARRIS = 1
    BILATERAL = 2
    LAPLACIAN = 3

    apps_list = set([])
    apps_list.add(NONE)
    apps_list.add(UNSHARP)
    apps_list.add(HARRIS)
    apps_list.add(BILATERAL)
    apps_list.add(LAPLACIAN)

    name = {}
    name[NONE] = ''
    name[UNSHARP] = 'Unsharp Mask'
    name[HARRIS] = 'Harris Corner'
    name[BILATERAL] = 'Bilateral Grid'
    name[LAPLACIAN] = 'Local Laplacian'

    file_name = {}
    file_name[NONE] = ''
    file_name[UNSHARP] = 'unsharp'
    file_name[HARRIS] = 'harris'
    file_name[BILATERAL] = 'bilateral'
    file_name[LAPLACIAN] = 'laplacian'

class ModeType():
    QUIT = -3
    PREV = -2
    NONE = -1
    CURRENT = 0
    CV2 = 1
    P_NAIVE = 2
    P_OPT = 3
    NUMBA = 4
    PIL = 5

    modes_list = set([])
    modes_list.add(CV2)
    modes_list.add(P_NAIVE)
    modes_list.add(P_OPT)
    modes_list.add(NUMBA)
    modes_list.add(PIL)

    all_modes = set([])
    all_modes.add(QUIT)
    all_modes.add(PREV)
    all_modes.add(NONE)
    all_modes.add(CURRENT)
    all_modes.union(modes_list)

    name = {}
    name[QUIT] = ''
    name[NONE] = ''
    name[CV2] = 'OpenCV'
    name[P_NAIVE] = 'PolyMage Naive'
    name[P_OPT] = 'PolyMage Opt'
    name[NUMBA] = 'Numpy + Numba'
    name[PIL] = 'Pillow (PIL)'

    @staticmethod
    def get_file_suffix(mode_id):
        if mode_id == ModeType.P_NAIVE:
            return '_naive'
        elif mode_id == ModeType.P_OPT:
            return ''

    @staticmethod
    def get_file_name(app_file, mode_id):
        suffix = ModeType.get_file_suffix(mode_id)
        file_name = app_file + suffix
        return file_name

    @staticmethod
    def get_func_name(lib, mode_id):
        if mode_id == ModeType.P_NAIVE:
            return lib.pipeline_naive
        elif mode_id == ModeType.P_OPT:
            return lib.pipeline_opt


class Mode:
    def __init__(self, _mode_id, _app_file, _py_func, _has_lib=False):
        modes_list = ModeType.modes_list.union([ModeType.NONE])
        assert _mode_id in modes_list, \
                "Mode type must be one of: %s" % modes_list + \
                "\nGiven: %s" % _mode_id
        self._mode_id = _mode_id
        self._app_file = _app_file
        self._py_func = _py_func
        self._has_lib = _has_lib

        self._lib_file = None
        self._lib = None
        self._lib_func = None

        self._init_libs()
        self._init_timers_frames()

    @property
    def mode_id(self):
        return self._mode_id
    @property
    def app_file(self):
        return self._app_file
    @property
    def py_func(self):
        return self._py_func
    @property
    def has_lib(self):
        return self._has_lib

    @property
    def lib_file(self):
        return self._lib_file

    @property
    def frames(self):
        return self._frames
    @property
    def time_spent(self):
        return self._total_time_spent
    @property
    def process_time(self):
        return self._process_time

    def _init_libs(self):
        self._lib_file = None
        self._lib = None
        self._lib_func = None
        if self.has_lib:
            # load the corresponding shared library
            file_name = ModeType.get_file_name(self.app_file, self.mode_id)
            self._lib_file = file_name+".so"
            mode_lib = ctypes.cdll.LoadLibrary(self._lib_file)
            # init memory pool pool
            mode_lib.pool_init()
            # set appropriate function names
            self._lib_func = ModeType.get_func_name(mode_lib, self.mode_id)
            self._lib = mode_lib
        return

    def _init_timers_frames(self):
        self._frames = 0
        self._process_time = 0.0
        self._total_time_spent = 0.0
        self._timeron = False
        return

    def _start_clock(self):
        if self._timeron == False:
            self._start_time = clock()
            self._timeron = True
        return

    def _stop_clock(self):
        if self._timeron == True:
            self._end_time = clock()
            self._process_time = self._end_time*1000 - self._start_time*1000
            self._timeron = False
        return

    def _frames_update(self):
        self._frames += 1
        return
    def _time_update(self, _time):
        self._total_time_spent += _time
        return

    def process_frame(self, frame):
        self._start_clock()
        res = self.py_func(frame, self._lib_func)
        self._stop_clock()

        self._frames_update()
        self._time_update(self._process_time)
        return res

    def avg_time_spent(self):
        return self.time_spent / self.frames

    def destroy(self):
        if self.has_lib:
            self._lib.pool_destroy()


class App:
    def __init__(self, _app_id, _dir, _mode_ids, _lib_mode_ids, _py_func_map):
        assert _app_id in AppType.apps_list, \
                "App ID must be one of %s" % AppType.apps_list
        self._app_id = _app_id
        self._name = AppType.name[_app_id]

        if self._app_id != AppType.NONE:
            assert os.path.exists(_dir), "%s : path does not exist" % _dir
            assert (len(_mode_ids) > 0), \
                "Apps can be run in at least one recognized mode"

        self._dir = _dir
        self._file = str(self._dir) + "/" + AppType.file_name[_app_id]

        assert self._file not in ["/"], \
            "path %s sounds fishy -_o" % self._file

        if _lib_mode_ids == None:
            _lib_mode_ids = []
        assert set(_lib_mode_ids).issubset(set(_mode_ids)), \
            "_lib_mode_ids (%s) should be a subset of " + \
            "allowed _mode_ids (%s)" % _lib_mode_ids % _mode_ids

        self._init_modes(_mode_ids, _lib_mode_ids, _py_func_map)
        self._set_initial_mode()

        self._process_time = 0.0

    @property
    def id_(self):
        return self._app_id
    @property
    def name(self):
        return self._name
    @property
    def file_name(self):
        return self._file
    @property
    def modes(self):
        return self._modes
    @property
    def current_mode(self):
        return self._cur_mode_id
    @property
    def previous_mode(self):
        return self._prev_mode_id
    @property
    def process_time(self):
        return self._process_time

    def _init_modes(self, _mode_ids, _lib_mode_ids, _py_func_map):
        self._modes = {}
        for mode_id in _mode_ids:
            has_lib = mode_id in _lib_mode_ids
            assert mode_id in _py_func_map, \
                "No frame processor function specified for mode %s" % mode_id
            py_func = _py_func_map[mode_id]
            self._modes[mode_id] = Mode(mode_id, self.file_name, py_func, has_lib)
        return

    def _set_initial_mode(self):
        mode_id = self._modes.keys()[0]
        self._cur_mode_id = mode_id
        self._prev_mode_id = mode_id
        return

    def process_frame(self, frame):
        if self.id_ == AppType.NONE:
            result_frame = frame
        else:
            cur_mode = self.modes[self._cur_mode_id]
            result_frame = cur_mode.process_frame(frame)
            self._process_time = cur_mode.process_time
        return result_frame

    def switch_mode(self, mode_id):
        if self.id_ == AppType.NONE:
            return
        if mode_id == ModeType.PREV:
            print "Go To PREV"
            temp = self._cur_mode_id
            self._cur_mode_id = self._prev_mode_id
            self._prev_mode_id = temp
        elif mode_id in self._modes:
            self._prev_mode_id = self._cur_mode_id
            self._cur_mode_id = mode_id

        print "----"
        print "cur:", self._cur_mode_id
        print "prev:", self._prev_mode_id
        return

    def destroy(self):
        for mode_id in self.modes:
            self.modes[mode_id].destroy()
        return

class VideoProcessor:
    def __init__(self, _apps_map, _file_name, _key_bind):
        self._apps_map = _apps_map
        self._file_name = _file_name
        self._key_bind = _key_bind

        self.set_initial_app(AppType.NONE)

        self._init_video()

    @property
    def apps_map(self):
        return self._apps_map
    @property
    def file_name(self):
        return self._file_name
    @property
    def key_bind(self):
        return self._key_bind
    @property
    def current_app(self):
        return self._current_app

    @property
    def frame_height(self):
        return self._frame_height
    @property
    def frame_width(self):
        return self._frame_width
    @property
    def frame_channels(self):
        return self._frame_channels

    @property
    def timer(self):
        return self._timer

    @property
    def process_time(self):
        return self._process_time

    def _switch_app(self, app_id):
        if app_id in self.apps_map:
            self._current_app = app_id
        return

    def set_initial_app(self, app_id):
        self._current_app = AppType.NONE
        self._switch_app(app_id)
        return

    def _init_video(self):
        # get video stream capturer
        assert os.path.exists(self.file_name), \
                "File (%s) not found" % self.file_name
        self._cap = VideoCapture(self.file_name)

        return

    def _update_frame(self):
        ret, self._cur_frame = self._cap.read()
        return

    def _process_frame(self, frame):
        app = self.apps_map[self.current_app]
        result_frame = app.process_frame(frame)
        result_frame = self._gen_output_frame(result_frame)
        return result_frame

    def _get_key_stroke(self):
        ch = 0xFF & waitKey(1)
        return ch

    def _get_next_change(self):
        '''
        Get the user's keystroke and update the changes in app and/or mode.
        Return bool if the next mode is to quit the video processing
        '''
        ch = self._get_key_stroke()
        if ch in self.key_bind[0]:
            app_id = self.current_app
            app = self.apps_map[app_id]
            mode_id = app.current_mode

            next_change = self.key_bind[((app_id, mode_id), ch)]
        else:
            return False

        # if app change is asked
        if next_change[0] != app_id:
            self._switch_app(next_change[0])
            return False
        else:
            # if app has not changed, look for the next mode change
            next_mode = next_change[1]
            if next_mode not in [ModeType.CURRENT, ModeType.QUIT]:
                app.switch_mode(next_mode)
            return next_mode == ModeType.QUIT

    def process(self, max_frames=0):
        stop = False
        if max_frames <= 0:
            # no upper limit
            while(self._cap.isOpened() and not stop):
                self._update_frame()
                result_frame = self._process_frame(self._cur_frame)
                self._display_output(result_frame)
                stop = self._get_next_change()
        else:
            # min of video #frames and #max_frames
            while(self._cap.isOpened() and \
                    self.frames < max_frames and not stop):
                self._update_frame()
                result_frame = self._process_frame(self._cur_frame)
                self._display_output(result_frame)
                stop = self._get_next_change()
        return

    def _gen_output_frame(self, frame):
        app = self.apps_map[self.current_app]

        label_start = (0, 0)
        label_end = (750, 150)
        colour = (255, 255, 255)
        rectangle(frame, label_start, label_end, colour, thickness=cv.CV_FILLED)

        # frame process time
        text1_start = (40, 40)
        text1 = "frame interval :  %.1f ms" % app.process_time
        draw_str(frame, text1_start, text1)

        # execution mode
        text2_start = (40, 80)
        text2 = "Pipeline       :  %s" % ModeType.name[app.current_mode]
        draw_str(frame, text2_start, text2)

        # app
        text3_start = (40, 120)
        text3 = "Benchmark      :  %s" % app.name
        draw_str(frame, text3_start, text3)

        return frame

    def _display_output(self, frame):
        frame = self._gen_output_frame(frame)
        imshow('Video', frame)
        return

    def report_stats(self):
        print "Average frame delays :"
        for app_id in self.apps_map:
            app = self.apps_map[app_id]
            for mode_id in app.modes:
                mode = app.modes[mode_id]
                if mode.frames:
                    avg = mode.avg_time_spent()
                    print "%s [%s] : %s ms" \
                        % (app.name, ModeType.name[mode_id], avg)
        return

    def finish(self):
        self._cap.release()
        destroyAllWindows()
        self.report_stats()
        pass
