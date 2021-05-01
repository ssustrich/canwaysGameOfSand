#![deny(clippy::all)]
#![forbid(unsafe_code)]

use env_logger::Env;
use log::{debug, error, info, warn};
use pixels::{Error, Pixels, SurfaceTexture};
use rand::Rng;
use winit::dpi::{LogicalPosition, LogicalSize, PhysicalSize};
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit_input_helper::WinitInputHelper;

const SCREEN_WIDTH: u32 = 10;
const SCREEN_HEIGHT: u32 = 10;
const INITIAL_FILL: f32 = 0.5;
const PARTICLETYPES: [&str; 3] = ["NONE", "SAND", "WATER"];
const GRAVITY: f32 = 9.8;
enum OBJECT_TYPES {
    EMPTY,
    WALL,
    SAND,
    WATER,
}
fn main() -> Result<(), Error> {
    env_logger::Builder::from_env(
        Env::default().default_filter_or("error,conways_gos=info"),
    )
    .init();
    log::info!("World dimensions: {}x{}", SCREEN_WIDTH, SCREEN_HEIGHT);

    let event_loop = EventLoop::new();
    let mut input = WinitInputHelper::new();
    let (window, p_width, p_height, mut _hidpi_factor) =
        create_window("Conway's Game of Sand", &event_loop);

    let surface_texture = SurfaceTexture::new(p_width, p_height, &window);

    //let mut frame = ConwayGrid::new_random(SCREEN_WIDTH as usize, SCREEN_HEIGHT as usize);
    let mut frame = ConwayGrid::new_empty(SCREEN_WIDTH as usize, SCREEN_HEIGHT as usize);
    let mut pixels = Pixels::new(SCREEN_WIDTH, SCREEN_HEIGHT, surface_texture)?;
    let mut paused = false;

    let mut draw_state: Option<bool> = None;

    event_loop.run(move |event, _, control_flow| {
        // The one and only event that winit_input_helper doesn't have for us...
        if let Event::RedrawRequested(_) = event {
            frame.draw(pixels.get_frame());
            if pixels
                .render()
                .map_err(|e| error!("pixels.render() failed: {}", e))
                .is_err()
            {
                *control_flow = ControlFlow::Exit;
                return;
            }
        }

        // For everything else, for let winit_input_helper collect events to build its state.
        // It returns `true` when it is time to update our game state and request a redraw.
        if input.update(&event) {
            // Close events
            if input.key_pressed(VirtualKeyCode::Escape) || input.quit() {
                *control_flow = ControlFlow::Exit;
                return;
            }
            if input.key_pressed(VirtualKeyCode::P) {
                paused = !paused;
            }
            if input.key_pressed(VirtualKeyCode::Space) {
                // Space is frame-step, so ensure we're paused
                paused = true;
            }
            if input.key_pressed(VirtualKeyCode::R) {
                frame.randomize();
            }
            // Handle mouse. This is a bit involved since support some simple
            // line drawing (mostly because it makes nice looking patterns).
            let (mouse_cell, mouse_prev_cell) = input
                .mouse()
                .map(|(mx, my)| {
                    let (dx, dy) = input.mouse_diff();
                    let prev_x = mx - dx;
                    let prev_y = my - dy;

                    let (mx_i, my_i) = pixels
                        .window_pos_to_pixel((mx, my))
                        .unwrap_or_else(|pos| pixels.clamp_pixel_pos(pos));

                    let (px_i, py_i) = pixels
                        .window_pos_to_pixel((prev_x, prev_y))
                        .unwrap_or_else(|pos| pixels.clamp_pixel_pos(pos));

                    (
                        (mx_i as isize, my_i as isize),
                        (px_i as isize, py_i as isize),
                    )
                })
                .unwrap_or_default();

            if input.mouse_pressed(0) {
                debug!("Mouse click at {:?}", mouse_cell);
                draw_state = Some(frame.toggle(mouse_cell.0, mouse_cell.1));
            } else if let Some(draw_alive) = draw_state {
                let release = input.mouse_released(0);
                let held = input.mouse_held(0);
                debug!("Draw at {:?} => {:?}", mouse_prev_cell, mouse_cell);
                debug!("Mouse held {:?}, release {:?}", held, release);
                // If they either released (finishing the drawing) or are still
                // in the middle of drawing, keep going.
                if release || held {
                    debug!("Draw line of {:?}", draw_alive);
                    frame.set_line(
                        mouse_prev_cell.0,
                        mouse_prev_cell.1,
                        mouse_cell.0,
                        mouse_cell.1,
                        draw_alive,
                    );
                }
                // If they let go or are otherwise not clicking anymore, stop drawing.
                if release || !held {
                    debug!("Draw end");
                    draw_state = None;
                }
            }
            // Adjust high DPI factor
            if let Some(factor) = input.scale_factor_changed() {
                _hidpi_factor = factor;
            }
            // Resize the window
            if let Some(size) = input.window_resized() {
                pixels.resize(size.width, size.height);
            }
            if !paused || input.key_pressed(VirtualKeyCode::Space) {
                frame.update();
            }
            window.request_redraw();
        }
    });
}

// COPYPASTE: ideally this could be shared.

//Methods for managaing meta game processes. Like opening a window or taking user input etc
fn create_window(
    title: &str,
    event_loop: &EventLoop<()>,
) -> (winit::window::Window, u32, u32, f64) {
    // Create a hidden window so we can estimate a good default window size
    let window = winit::window::WindowBuilder::new()
        .with_visible(false)
        .with_title(title)
        .build(&event_loop)
        .unwrap();
    let hidpi_factor = window.scale_factor();

    // Get dimensions
    let width = SCREEN_WIDTH as f64;
    let height = SCREEN_HEIGHT as f64;
    let (monitor_width, monitor_height) = {
        if let Some(monitor) = window.current_monitor() {
            let size = monitor.size().to_logical(hidpi_factor);
            (size.width, size.height)
        } else {
            (width, height)
        }
    };
    let scale = (monitor_height / height * 2.0 / 3.0).round().max(1.0);

    // Resize, center, and display the window
    let min_size: winit::dpi::LogicalSize<f64> =
        PhysicalSize::new(width, height).to_logical(hidpi_factor);
    let default_size = LogicalSize::new(width * scale, height * scale);
    let center = LogicalPosition::new(
        (monitor_width - width * scale) / 2.0,
        (monitor_height - height * scale) / 2.0,
    );
    window.set_inner_size(default_size);
    window.set_min_inner_size(Some(min_size));
    window.set_outer_position(center);
    window.set_visible(true);

    let size = default_size.to_physical::<f64>(hidpi_factor);

    (
        window,
        size.width.round() as u32,
        size.height.round() as u32,
        hidpi_factor,
    )
}

/// Generate a pseudorandom seed for the game's PRNG.
fn generate_seed() -> (u64, u64) {
    use byteorder::{ByteOrder, NativeEndian};
    use getrandom::getrandom;

    let mut seed = [0_u8; 16];

    getrandom(&mut seed).expect("failed to getrandom");

    (
        NativeEndian::read_u64(&seed[0..8]),
        NativeEndian::read_u64(&seed[8..16]),
    )
}

#[derive(Clone, Copy, Debug, Default)]
/// The most basic element in teh game
struct Particle {
    p_type: usize,
    active: bool,
    already_updated: bool,
    velocity: f32,
}

impl Particle {
    fn new(p_type: usize, active: bool) -> Self {
        Self {
            p_type,
            active,
            already_updated: false,
            velocity: 1.0,
        }
    }

    #[must_use]
    fn next_state(mut self, active: bool) -> Self {
        self.active = active;
        self
    }

    fn set_active(&mut self, active: bool) {
        *self = self.next_state(active);
    }
}

#[derive(Clone, Debug)]
struct ConwayGrid {
    particles: Vec<Particle>,
    width: usize,
    height: usize,
    // Should always be the same size as `cells`. When updating, we read from
    // `cells` and write to `scratch_cells`, then swap. Otherwise it's not in
    // use, and `cells` should be updated directly.
    scratch_cells: Vec<Particle>,
}

impl ConwayGrid {
    fn new_empty(width: usize, height: usize) -> Self {
        assert!(width != 0 && height != 0);
        let size = width.checked_mul(height).expect("too big");
        Self {
            particles: vec![Particle::default(); size],
            scratch_cells: vec![Particle::default(); size],
            width,
            height,
        }
    }

    fn new_random(width: usize, height: usize) -> Self {
        let mut result = Self::new_empty(width, height);
        result.randomize();
        result
    }

    fn randomize(&mut self) {
        let mut rng = rand::thread_rng();
        for x in 0..10 {
            let x = rng.gen_range(0..(SCREEN_HEIGHT * SCREEN_HEIGHT - 10)) as usize;
            self.particles[x].set_active(true);
            self.particles[x].p_type = 1;
            self.particles[x].already_updated = false;
        }

        // let mut rng: randomize::PCG32 = generate_seed().into();
        // for c in self.particles.iter_mut() {
        //     let alive = randomize::f32_half_open_right(rng.next_u32()) > INITIAL_FILL;
        //     let mut rng1 = rand::thread_rng();
        //     //println!("Integer: {}", rng1.gen_range(0..PARTICLETYPES.len());
        //     *c = Particle::new(rng1.gen_range(0..PARTICLETYPES.len()), true);
        // }
        // run a few simulation iterations for aesthetics (If we don't, the
        // noise is ugly)
        for _ in 0..3 {
            self.update();
        }
    }

    fn update(&mut self) {
        let mut counter = 0;
        for idx in (0..self.particles.len()).rev() {
            //let neibs = self.count_neibs(x, y);
            //println!("Checking for alive cell at index {:?}", self.particles[idx]);
            //self.printCrazy8(v, idx);

            if !self.particles[idx].already_updated && self.particles[idx].active && self.particles[idx].p_type ==1  {
                
                //if at base level set to inactive go to next particle
                //if there is not a particle below it then move it down one akd keep it active
                //if there a dir particle below it can it go left or right
                //if so then go bl or br and set to inactince
                //else stop and set to inactive

                
                counter+=1;
                println!("counter is now {} ", counter);
                //println!("executing loop counter {}",counter);
                log::debug!("{:?}", self.particles[idx]);
                //check to see if we can move down
                let mut v: Vec<isize> = self.getEightNeighbors(idx);
                let mut bi = v[2] as isize;
                //we hit the bottom
                if bi == -1 {
                    self.particles[idx].active = false;
                    self.particles[idx].already_updated = true;
                } else if  self.particles[bi as usize].p_type == 0{
                    self.particles[idx].active = false;
                    self.particles[idx].already_updated = true;
                    self.particles[idx].p_type =0;
                    self.particles[bi as usize].active = true;
                    self.particles[bi as usize].already_updated = false;
                    self.particles[bi as usize].p_type =1;
                } else if self.particles[bi as usize].p_type == 1{
                    let mut bl = v[3] as isize;
                    let mut br = v[1] as isize;
                    if rand::random() {
                            bl ^= br;
                            br ^= bl;
                            bl ^= br;
                    }
                    if  bl > -1 && self.particles[bl as usize].p_type == 0  {
                        self.particles[idx].active = false;
                        self.particles[idx].already_updated = true;
                        self.particles[idx].p_type =0;
                        self.particles[bl as usize].active = true;
                        self.particles[bl as usize].already_updated = false;
                        self.particles[bl as usize].p_type =1;
                    }
                    else if  br > -1 && self.particles[br as usize].p_type == 0{
                        self.particles[idx].active = false;
                        self.particles[idx].already_updated = true;
                        self.particles[idx].p_type =0;
                        self.particles[br as usize].active = true;
                        self.particles[br as usize].already_updated = false;
                        self.particles[br as usize].p_type =1;
                    }
                    else if  self.particles[bi as usize].active == true {
                        self.particles[idx].active = false;
                        self.particles[idx].already_updated = true;
                        self.particles[idx].p_type =0;
                        self.particles[bi as usize].active = true;
                        self.particles[bi as usize].already_updated = false;
                        self.particles[bi as usize].p_type =1;
                    }else if  self.particles[bi as usize].active == false {
                         self.particles[idx].active = false;
                        self.particles[idx].already_updated = true;
                    }
                }
                
            }    
                
                
                
                
                
                
                
                
                



        //             if self.particles[bi as usize].p_type == 0 {
        //                 self.particles[idx].already_updated = true;
        //                 self.particles[idx].active = false;
        //                 self.particles[idx].p_type = 0;
        //                 self.particles[bi as usize].active = true;
        //                 self.particles[bi as usize].p_type = 1;
        //             } 
        //             else if self.particles[bi as usize].p_type == 1 {
        //            //     self.particles[idx].active = false;
        //            //     self.particles[idx].already_updated = true;
        //
        //             // if rand::random() {
        //             //     bl ^= br;
        //             //     br ^= bl;
        //             //     bl ^= br;
        //             // }
        //
        //             // //check bl  (which may be swapped)
        //               if bl > -1 && self.particles[bl as usize].p_type == 0{
        //                   self.particles[idx].already_updated = true;
        //                   self.particles[idx].active = false;
        //                   self.particles[idx].p_type =0;
        //                   self.particles[bl as usize].already_updated = true;
        //                   self.particles[bl as usize].active = true;
        //                   self.particles[bl as usize].p_type = 1;
        //               }
        //             // else if br > -1 && self.particles[br as usize].p_type == 0{
        //             //     self.particles[idx].already_updated = true;
        //             //     self.particles[idx].active = false;
        //             //     self.particles[idx].p_type =0;
        //             //     self.particles[br as usize].already_updated = true;
        //             //     self.particles[br as usize].active = true;
        //             //     self.particles[br as usize].p_type = 1;
        //             // }
        //         }
        //              else{
        //              self.particles[idx].already_updated = true;
        //              self.particles[idx].active = false;
        //              }
        //         } 
        //         // else {
        //         //     self.particles[idx].already_updated = true;
        //         //     self.particles[idx].active = true;
        //         //     // self.particles[bi as usize].active = true;
        //         // }
        //     }
        //}
        // for y in 0..self.height {
        //     for x in 0..self.width {
        //         let idx = x + y * self.width;
        //         self.particles[idx].already_updated = false;
        //     }
        // }
        }
    }

    fn toggle(&mut self, x: isize, y: isize) -> bool {
        println!("Toggle called");
        if let Some(i) = self.grid_idx(x, y) {
            println!("On inded {}", i);
            let was_active = self.particles[i].active;
            println!("its active state is{}", was_active);
            self.particles[i].set_active(!was_active);
            println!("Now its active state is {}",  self.particles[i].active);

            if !was_active {
                self.particles[i].p_type = 1;
                self.particles[i].already_updated = !true;
                println!("p_type is {}",  self.particles[i].p_type);
            }
            else{
                self.particles[i].p_type = 0;
                self.particles[i].already_updated = true;
                println!("p type is  {}",  self.particles[i].p_type);
            }
            
            println!("We are returning {}", !was_active);

            !was_active
        } else {
            println!("Do we come here often?");
            false
        }
    }

    fn draw(&self, screen: &mut [u8]) {
        debug_assert_eq!(screen.len(), 4 * self.particles.len());
        for (c, pix) in self.particles.iter().zip(screen.chunks_exact_mut(4)) {
            let mut color = [0, 0x00, 0x00, 0x00];
            if c.p_type == 1 {
                color = [0, 0xff, 0xff, 0xff];
            }
            if c.p_type == 0 {
                color = [0, 0x00, 0x00, 0x00];
            }

            // let color = if c.active {
            //     [0, 0xff, 0xff, 0xff]
            // } else {
            //     [0, 0, 0, 0xff]
            // };
            pix.copy_from_slice(&color);
        }
    }

    fn set_line(&mut self, x0: isize, y0: isize, x1: isize, y1: isize, active: bool) {
        // probably should do sutherland-hodgeman if this were more serious.
        // instead just clamp the start pos, and draw until moving towards the
        // end pos takes us out of bounds.
        let x0 = x0.max(0).min(self.width as isize);
        let y0 = y0.max(0).min(self.height as isize);
        for (x, y) in line_drawing::Bresenham::new((x0, y0), (x1, y1)) {
            if let Some(i) = self.grid_idx(x, y) {
                self.particles[i].set_active(active);
                self.particles[i].p_type = 1;
            } else {
                break;
            }
        }
    }

    fn getXYfromInx(&self, idx: usize) -> (usize, usize) {
        let row: usize = idx % self.width;
        let column: usize = idx / self.width;
        (row, column)
    }

    /* Given an index in the array of X this function will return the index on the 8
    neighbors in an array of len 8 where array[0] is the cell to the immediate right of ?
    and the continue in a clowise fashion. If the cell is touching an edge of the game board the value
    for neighbors that are off the board is -1

    [5]  [6]  [7]

    [4]   X   [0]

    [3]  [2]  [1]


    Example 1) Cell X is in the middle of a 3x3 game board

    idx = 4 and the board would look like this
    0  1  2
    3  4  5
    6  7  8
    The returned vector is <5,8,7,6,3,0,1,2>

    Example 2) Cell X is in the upper left corner of a 3x3 board
    idx = 0 and the board would look like this, remember -1 indicataes a wall or edge
    -1 -1 -1
    -1  0  1
    -1  3  4
    The returned vector is <1,4,3,-1,-1,-1,-1,-1>

    */
    fn getEightNeighbors(&self, idx: usize) -> Vec<isize> {
        let coord = self.getXYfromInx(idx);
        let mut v: Vec<isize> = vec![-1; 8];
        log::info!(
            "[getEightNeighbors] Called with indx:{}, which maps to x,y:{:?}",
            idx,
            coord
        );
        //println!("{}{:?}",idx, coord );

        v[0] = if coord.0 == self.width - 1 {
            -1
        } else {
            (idx + 1) as isize
        };
        v[1] = if coord.0 == self.width - 1 || coord.1 == self.height - 1 {
            -1
        } else {
            (idx + 1 + self.width) as isize
        };
        v[2] = if coord.1 == self.height - 1 {
            -1
        } else {
            (idx + self.width) as isize
        };
        v[3] = if coord.0 == 0 || coord.1 == self.height - 1 {
            -1
        } else {
            (idx + self.width - 1) as isize
        };
        v[4] = if coord.0 == 0 { -1 } else { (idx - 1) as isize };
        v[5] = if coord.1 == 0 || coord.0 == 0 {
            -1
        } else {
            (idx - 1 - self.width) as isize
        };
        v[6] = if coord.1 == 0 {
            -1
        } else {
            (idx - self.width) as isize
        };
        v[7] = if coord.0 == self.width - 1 || coord.1 == 0 {
            -1
        } else {
            (idx + 1 - self.width) as isize
        };
        log::debug!("[getEightNeighbors] {:?}", v);
        v
    }

    fn printCrazy8(&self, x: Vec<isize>, idx: usize) {
        println!();
        println!("{} {} {}", x[5], x[6], x[7]);
        println!("{} {} {}", x[4], idx, x[0]);
        println!("{} {} {}", x[3], x[2], x[1]);
        println!();
    }

    fn grid_idx<I: std::convert::TryInto<usize>>(&self, x: I, y: I) -> Option<usize> {
        if let (Ok(x), Ok(y)) = (x.try_into(), y.try_into()) {
            if x < self.width && y < self.height {
                Some(x + y * self.width)
            } else {
                None
            }
        } else {
            None
        }
    }
}

/*
fn main() {

    let a:i32 = 0;
    let b:i32 = 100;
    let width:i32 = 10;
    let mut c:i8  = 0;
    for x in a..b{
        if x<width{
             c = (x.overflowing_sub(width)).0 as i8;
        }
        else{
            c = (x-width) as i8;
            }
        println!("{}", c);
    }
    }

    */
