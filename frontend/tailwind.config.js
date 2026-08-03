/** @type {import('tailwindcss').Config} */
export default {
  content: ['./index.html', './src/**/*.{js,ts,jsx,tsx}'],
  theme: {
    extend: {
      colors: {
        // Navy / steel / grey, matching the house style of the meeting decks
        // and the html_dphil explainers. No purple.
        ink: {
          900: '#0f1e2e',
          800: '#16293c',
          700: '#1e3a52',
          600: '#2b5070',
          500: '#3d6b91',
          400: '#5b89ae',
          300: '#8badc9',
          200: '#c3d6e4',
          100: '#e4edf4',
          50: '#f4f8fb',
        },
        // Element classes, fixed by the `re` column in annotated.tsv.
        element: {
          enhancer: '#c2703d',
          ctcf: '#4a7c59',
          promoter: '#2b5070',
        },
      },
      fontFamily: {
        sans: ['Inter', 'ui-sans-serif', 'system-ui', 'sans-serif'],
        mono: ['ui-monospace', 'SFMono-Regular', 'Menlo', 'monospace'],
      },
    },
  },
  plugins: [],
}
