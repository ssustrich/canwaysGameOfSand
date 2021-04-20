#![deny(clippy::all)]
#![forbid(unsafe_code)]

use log::{debug, error};
use pixels::{Error, Pixels, SurfaceTexture};
use rand::Rng;
use winit::dpi::{LogicalPosition, LogicalSize, PhysicalSize};
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit_input_helper::WinitInputHelper;

const SCREEN_WIDTH: u32 = 5;
const SCREEN_HEIGHT: u32 = 5;
const INITIAL_FILL: f32 = 0.9;
const PARTICLETYPES: [&str; 3] = ["NONE", "SAND", "WATER"];
const GRAVITY: f32 = 9.8;
enum OBJECT_TYPES {
    EMPTY,
    WALL,
    SAND,
    WATER,
}

fn main() -> Result<(), Error> {
    env_logger::init();
    let event_loop = EventLoop::new();
    let mut input = WinitInputHelper::new();
    let (window, p_width, p_height, mut _hidpi_factor) =
        create_window("Conway's Game of Sand", &event_loop);

    let surface_texture = SurfaceTexture::new(p_width, p_height, &window);

    let mut frame = ConwayGrid::new_random(SCREEN_WIDTH as usize, SCREEN_HEIGHT as usize);
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
                pixels.resize_surface(size.width, size.height);
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
struct Particle {
    p_type: usize,
    active: bool,
    already_updated: bool,
    heat: u8,
    velocity: f32,
}

impl Particle {
    fn new(p_type: usize, active: bool) -> Self {
        Self {
            p_type,
            active,
            already_updated: false,
            heat: 0,
            velocity: 0.0,
        }
    }

    #[must_use]
    fn next_state(mut self, active: bool) -> Self {
        self.active = active;
        if self.active {
            self.heat = 0;
        } else {
            self.heat = 0;
        }
        self
    }

    fn set_active(&mut self, active: bool) {
        *self = self.next_state(active);
    }

    fn cool_off(&mut self, decay: f32) {
        if !self.active {
            let heat = (self.heat as f32 * decay).min(255.0).max(0.0);
            assert!(heat.is_finite());
            self.heat = heat as u8;
        }
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
        let mut rng: randomize::PCG32 = generate_seed().into();
        for c in self.particles.iter_mut() {
            let alive = randomize::f32_half_open_right(rng.next_u32()) > INITIAL_FILL;
            let mut rng1 = rand::thread_rng();
            //println!("Integer: {}", rng1.gen_range(0..PARTICLETYPES.len());
            *c = Particle::new(rng1.gen_range(0..PARTICLETYPES.len()), alive);
        }
        // run a few simulation iterations for aesthetics (If we don't, the
        // noise is ugly)
        for _ in 0..3 {
            self.update();
        }
        // Smooth out noise in the heatmap that would remain for a while
        for c in self.particles.iter_mut() {
            c.cool_off(0.4);
        }
    }

    fn update(&mut self) {
        for y in 0..self.height {
            for x in 0..self.width {
                //let neibs = self.count_neibs(x, y);
                let idx = x + y * self.width;
                println!("Checking for alive cell [{}{}] at index {}", x, y, idx,);
                self.getEightNeighbors(idx);
                if self.particles[idx].active
                    && idx + self.width < self.particles.len()
                    && self.particles[idx + self.width].already_updated == false
                    && self.particles[idx].already_updated == false
                {
                    //  println!("Cell at index {} is alive", idx);
                    //    println!("Killing cell and dropping particle to the cell bellow us at {}", idx+self.width);
                    if self.particles[idx + self.width].active == false {
                        //let num = rand::thread_rng().gen_range(0..2);
                        self.particles[idx + self.width].active = true;
                        self.particles[idx + self.width].already_updated = true;
                        self.particles[idx].active = false;
                    }
                    // else if (idx + self.width -1) < self.cells.len()
                    // && (idx + self.width) % self.width != 0
                    // && self.cells[idx + self.width -1].alive == false
                    // {
                    //   self.cells[idx + self.width -1].alive = true;
                    //   self.cells[idx + self.width].already_updated = true;
                    //   self.cells[idx].alive = false;
                    //  }
                    //  else if (idx + self.width +1) < self.cells.len()
                    //  && (idx + self.width) % self.width != 0
                    //  && self.cells[idx + self.width +1].alive == false
                    //  {
                    //    self.cells[idx + self.width +1].alive = true;
                    //    self.cells[idx + self.width].already_updated = true;
                    //    self.cells[idx].alive = false;
                    //  }
                }
            }
        }
        //   std::mem::swap(&mut self.scratch_cells, &mut self.cells);
        //        println!("Do we get here?");
        for y in 0..self.height {
            for x in 0..self.width {
                let idx = x + y * self.width;
                self.particles[idx].already_updated = false;
            }
        }
    }

    fn toggle(&mut self, x: isize, y: isize) -> bool {
        if let Some(i) = self.grid_idx(x, y) {
            let was_alive = self.particles[i].active;
            self.particles[i].set_active(!was_alive);
            !was_alive
        } else {
            false
        }
    }

    fn draw(&self, screen: &mut [u8]) {
        debug_assert_eq!(screen.len(), 4 * self.particles.len());
        for (c, pix) in self.particles.iter().zip(screen.chunks_exact_mut(4)) {
            let color = if c.active {
                [0, 0xff, 0xff, 0xff]
            } else {
                [0, 0, c.heat, 0xff]
            };
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
            } else {
                break;
            }
        }
    }

    fn getXYfromInx(&self, idx: usize) -> (usize, usize) {
        let row: usize = idx / self.width;
        let column: usize = idx / self.width;
        (row, column)
    }

    fn getEightNeighbors(&self, idx: usize) -> Vec<isize> {
        if idx == 6 {
            println!("The width is {}", self.width);
        }
        println!("The width is {}", self.width);
        // let (x,y) = self.getXYfromInx(idx);
        let mut v: Vec<isize> = vec![-1; 8];
        // let mut x_to_check: usize = x as i32 - (-1);
        // let mut y_to_check: usize = y as i32 - (-1);
        // for y in 0..8{
        //     for x in 0..8{
        //         v[x as usize +y as usize] = self.grid_idx(x_to_check as i32, y_to_check as i32);
        //         x_to_check +=1;
        //     }

        // }
        // v[0]= idx.wrapping_sub(self.width - 1);

        for x in 0..8 {
            if x == 0 {
                if idx <= self.width  || idx % self.width == 0  {
                    v[0] = -1;
                } else {
                    v[0] = (idx - self.width - 1) as isize;
                }
            } else if x == 1 {
                if idx <= self.width {
                    v[1] =  -1;
                } else {
                    v[1] = (idx - self.width) as isize;
                }
            } else if x == 2   || idx + 1 % self.width == 0  {
                if idx <= self.width {
                    v[2] =  -1;
                } else {
                    v[2] = (idx - self.width + 1) as isize;
                }
            } else if x == 3 {
                if idx % self.width == 0  {
                    v[3] =  -1;
                } else {
                    v[3] = (idx - 1) as isize;
                }
            } else if x == 4 {
                if idx + 1 % self.width == 0 {
                    v[4] =  -1;
                } else {
                    v[4] = (idx + 1) as isize;
                }
            } else if x == 5 {
                if idx + 1 % self.width == 0 {
                    v[5] =  -1;
                } else {
                    v[5] = (idx + self.width - 1) as isize;
                }
            } else if x == 6 {
                v[6] = (idx + self.width) as isize;
            } else if x == 7 {
                if idx + 1 % self.width == 0 {
                    v[7] =  -1;
                } else {
                    v[7] = (idx + self.width + 1) as isize;
                }
            }
        }
        // if idx <= self.width {
        //     v[0] = -1;
        //     v[1] = -1;
        //     v[2] = -1;
        // } else {
        //     v[0] = (idx - self.width - 1) as isize;
        //     v[1] = (idx - self.width) as isize;
        //     v[2] = (idx - self.width + 1) as isize;
        // }
        // if idx + 1 % self.width == 0 {
        //     v[0] = -1;
        //     v[3] = -1;
        //     v[5] = -1;
        // } else {
        //     v[0] = (idx - self.width - 1) as isize;
        //     v[3] = (idx - 1) as isize;
        //     v[5] = (idx + self.width - 1) as isize;
        // }

        // if idx + 1 % self.width == 0 {
        //     v[2] = (idx - self.width + 1) as isize;
        //     v[4] = (idx + 1) as isize;
        //     v[7] = (idx + self.width + 1) as isize;
        // } else {
        //     v[2] = (idx - self.width + 1) as isize;
        //     v[4] = (idx + 1) as isize;
        //     v[7] = (idx + self.width + 1) as isize;
        // }

        // v[6] = (idx + self.width) as isize;
        println!("The eight neighbors of {} are {:?}", idx, v);
        v
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
